% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RPRRRR10_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_jacobig_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:19:50
% EndTime: 2019-02-26 21:19:50
% DurationCPUTime: 0.06s
% Computational Cost: add. (25->13), mult. (77->31), div. (0->0), fcn. (111->10), ass. (0->22)
t140 = sin(pkin(6));
t145 = sin(qJ(1));
t153 = t140 * t145;
t147 = cos(qJ(1));
t152 = t140 * t147;
t138 = sin(pkin(13));
t151 = t145 * t138;
t141 = cos(pkin(13));
t150 = t145 * t141;
t149 = t147 * t138;
t148 = t147 * t141;
t146 = cos(qJ(3));
t144 = sin(qJ(3));
t143 = cos(pkin(6));
t142 = cos(pkin(7));
t139 = sin(pkin(7));
t137 = -t143 * t150 - t149;
t136 = t143 * t148 - t151;
t135 = -t143 * t139 * t146 + (-t141 * t142 * t146 + t138 * t144) * t140;
t134 = (-t143 * t151 + t148) * t144 + (-t137 * t142 - t139 * t153) * t146;
t133 = (t143 * t149 + t150) * t144 + (-t136 * t142 + t139 * t152) * t146;
t1 = [0, 0, -t137 * t139 + t142 * t153, t134, t134, 0; 0, 0, -t136 * t139 - t142 * t152, t133, t133, 0; 1, 0, -t139 * t140 * t141 + t142 * t143, t135, t135, 0;];
Jg_rot  = t1;
