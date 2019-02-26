% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PPPRRR1
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:38
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PPPRRR1_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_jacobig_rot_5_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:38:44
% EndTime: 2019-02-26 19:38:44
% DurationCPUTime: 0.07s
% Computational Cost: add. (49->26), mult. (144->58), div. (0->0), fcn. (199->14), ass. (0->33)
t133 = sin(pkin(12));
t142 = cos(pkin(6));
t152 = t133 * t142;
t135 = sin(pkin(7));
t136 = sin(pkin(6));
t151 = t135 * t136;
t150 = t135 * t142;
t141 = cos(pkin(7));
t149 = t136 * t141;
t138 = cos(pkin(13));
t148 = t138 * t141;
t139 = cos(pkin(12));
t147 = t139 * t142;
t132 = sin(pkin(13));
t127 = -t133 * t132 + t138 * t147;
t146 = t127 * t141 - t139 * t151;
t129 = -t139 * t132 - t138 * t152;
t145 = t129 * t141 + t133 * t151;
t144 = cos(qJ(4));
t143 = sin(qJ(4));
t140 = cos(pkin(8));
t137 = cos(pkin(14));
t134 = sin(pkin(8));
t131 = sin(pkin(14));
t130 = -t132 * t152 + t139 * t138;
t128 = t132 * t147 + t133 * t138;
t126 = -t138 * t151 + t142 * t141;
t125 = -t129 * t135 + t133 * t149;
t124 = -t127 * t135 - t139 * t149;
t123 = t137 * t150 + (-t131 * t132 + t137 * t148) * t136;
t122 = -t130 * t131 + t145 * t137;
t121 = -t128 * t131 + t146 * t137;
t1 = [0, 0, 0, -t122 * t134 + t125 * t140 (t130 * t137 + t145 * t131) * t143 + (-t122 * t140 - t125 * t134) * t144, 0; 0, 0, 0, -t121 * t134 + t124 * t140 (t128 * t137 + t146 * t131) * t143 + (-t121 * t140 - t124 * t134) * t144, 0; 0, 0, 0, -t123 * t134 + t126 * t140 (t136 * t132 * t137 + (t136 * t148 + t150) * t131) * t143 + (-t123 * t140 - t126 * t134) * t144, 0;];
Jg_rot  = t1;
