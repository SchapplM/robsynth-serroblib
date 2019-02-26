% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RPRPRR9_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_jacobig_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:53:24
% EndTime: 2019-02-26 20:53:24
% DurationCPUTime: 0.04s
% Computational Cost: add. (24->15), mult. (70->33), div. (0->0), fcn. (100->12), ass. (0->25)
t137 = sin(pkin(6));
t143 = sin(qJ(1));
t152 = t137 * t143;
t145 = cos(qJ(1));
t151 = t137 * t145;
t135 = sin(pkin(12));
t150 = t143 * t135;
t139 = cos(pkin(12));
t149 = t143 * t139;
t148 = t145 * t135;
t147 = t145 * t139;
t134 = sin(pkin(13));
t138 = cos(pkin(13));
t142 = sin(qJ(3));
t144 = cos(qJ(3));
t146 = -t134 * t142 + t138 * t144;
t141 = cos(pkin(6));
t140 = cos(pkin(7));
t136 = sin(pkin(7));
t133 = -t144 * t134 - t142 * t138;
t132 = -t141 * t149 - t148;
t131 = t141 * t147 - t150;
t130 = t146 * t140;
t129 = t146 * t136;
t1 = [0, 0, -t132 * t136 + t140 * t152, 0 -(-t141 * t150 + t147) * t133 - t132 * t130 - t129 * t152, 0; 0, 0, -t131 * t136 - t140 * t151, 0 -(t141 * t148 + t149) * t133 - t131 * t130 + t129 * t151, 0; 1, 0, -t137 * t139 * t136 + t141 * t140, 0, -t141 * t129 + (-t130 * t139 - t133 * t135) * t137, 0;];
Jg_rot  = t1;
