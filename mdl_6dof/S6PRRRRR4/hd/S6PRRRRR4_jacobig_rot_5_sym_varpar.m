% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRRR4
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRRRR4_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR4_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR4_jacobig_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:20:26
% EndTime: 2019-02-26 20:20:26
% DurationCPUTime: 0.05s
% Computational Cost: add. (26->14), mult. (79->33), div. (0->0), fcn. (116->10), ass. (0->21)
t136 = sin(pkin(13));
t138 = sin(pkin(6));
t150 = t136 * t138;
t137 = sin(pkin(7));
t149 = t137 * t138;
t139 = cos(pkin(13));
t148 = t139 * t138;
t141 = cos(pkin(6));
t143 = sin(qJ(2));
t147 = t141 * t143;
t145 = cos(qJ(2));
t146 = t141 * t145;
t144 = cos(qJ(3));
t142 = sin(qJ(3));
t140 = cos(pkin(7));
t135 = -t136 * t146 - t139 * t143;
t134 = -t136 * t143 + t139 * t146;
t133 = -t141 * t137 * t144 + (-t140 * t144 * t145 + t142 * t143) * t138;
t132 = (-t136 * t147 + t139 * t145) * t142 + (-t135 * t140 - t136 * t149) * t144;
t131 = (t136 * t145 + t139 * t147) * t142 + (-t134 * t140 + t137 * t148) * t144;
t1 = [0, t150, -t135 * t137 + t140 * t150, t132, t132, 0; 0, -t148, -t134 * t137 - t140 * t148, t131, t131, 0; 0, t141, t141 * t140 - t145 * t149, t133, t133, 0;];
Jg_rot  = t1;
