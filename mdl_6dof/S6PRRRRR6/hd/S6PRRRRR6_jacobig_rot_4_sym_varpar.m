% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRRRRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRRRR6_jacobig_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_jacobig_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_jacobig_rot_4_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:21:51
% EndTime: 2019-02-26 20:21:52
% DurationCPUTime: 0.05s
% Computational Cost: add. (24->17), mult. (69->39), div. (0->0), fcn. (101->12), ass. (0->23)
t142 = sin(pkin(14));
t145 = sin(pkin(6));
t158 = t142 * t145;
t144 = sin(pkin(7));
t157 = t144 * t145;
t146 = cos(pkin(14));
t156 = t146 * t145;
t149 = cos(pkin(6));
t151 = sin(qJ(2));
t155 = t149 * t151;
t153 = cos(qJ(2));
t154 = t149 * t153;
t152 = cos(qJ(3));
t150 = sin(qJ(3));
t148 = cos(pkin(7));
t147 = cos(pkin(8));
t143 = sin(pkin(8));
t141 = -t142 * t154 - t146 * t151;
t140 = -t142 * t151 + t146 * t154;
t139 = t149 * t148 - t153 * t157;
t138 = -t141 * t144 + t148 * t158;
t137 = -t140 * t144 - t148 * t156;
t1 = [0, t158, t138 -(-(-t142 * t155 + t146 * t153) * t150 + (t141 * t148 + t142 * t157) * t152) * t143 + t138 * t147, 0, 0; 0, -t156, t137 -(-(t142 * t153 + t146 * t155) * t150 + (t140 * t148 - t144 * t156) * t152) * t143 + t137 * t147, 0, 0; 0, t149, t139 -(t149 * t144 * t152 + (t148 * t152 * t153 - t150 * t151) * t145) * t143 + t139 * t147, 0, 0;];
Jg_rot  = t1;
