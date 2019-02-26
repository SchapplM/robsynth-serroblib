% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPRRR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:56
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRPRRR14_jacobig_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobig_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobig_rot_4_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:55:45
% EndTime: 2019-02-26 22:55:45
% DurationCPUTime: 0.05s
% Computational Cost: add. (19->17), mult. (54->37), div. (0->0), fcn. (78->12), ass. (0->21)
t150 = sin(pkin(6));
t156 = sin(qJ(1));
t164 = t156 * t150;
t155 = sin(qJ(2));
t163 = t156 * t155;
t157 = cos(qJ(2));
t162 = t156 * t157;
t158 = cos(qJ(1));
t161 = t158 * t150;
t160 = t158 * t155;
t159 = t158 * t157;
t154 = cos(pkin(6));
t153 = cos(pkin(7));
t152 = cos(pkin(8));
t151 = cos(pkin(14));
t149 = sin(pkin(7));
t148 = sin(pkin(8));
t147 = sin(pkin(14));
t146 = -t154 * t162 - t160;
t145 = t154 * t159 - t163;
t1 = [0, t164, 0 -(-(-t154 * t163 + t159) * t147 + (t146 * t153 + t149 * t164) * t151) * t148 + (-t146 * t149 + t153 * t164) * t152, 0, 0; 0, -t161, 0 -(-(t154 * t160 + t162) * t147 + (t145 * t153 - t149 * t161) * t151) * t148 + (-t145 * t149 - t153 * t161) * t152, 0, 0; 1, t154, 0 -(t154 * t149 * t151 + (t151 * t153 * t157 - t147 * t155) * t150) * t148 + (-t150 * t157 * t149 + t154 * t153) * t152, 0, 0;];
Jg_rot  = t1;
