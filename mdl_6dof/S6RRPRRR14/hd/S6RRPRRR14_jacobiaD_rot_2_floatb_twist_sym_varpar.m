% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 2 (0=Basis) von
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-12-10 18:38
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRR14_jacobiaD_rot_2_floatb_twist_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobiaD_rot_2_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14_jacobiaD_rot_2_floatb_twist_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobiaD_rot_2_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_2_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-10 18:38:20
% EndTime: 2018-12-10 18:38:21
% DurationCPUTime: 0.32s
% Computational Cost: add. (564->49), mult. (1035->121), div. (133->12), fcn. (1112->13), ass. (0->61)
t140 = pkin(6) + qJ(2);
t141 = pkin(6) - qJ(2);
t103 = cos(t141) / 0.2e1 + cos(t140) / 0.2e1;
t123 = sin(qJ(2));
t124 = sin(qJ(1));
t126 = cos(qJ(1));
t94 = t124 * t103 + t126 * t123;
t88 = t94 ^ 2;
t102 = sin(t140) / 0.2e1 - sin(t141) / 0.2e1;
t125 = cos(qJ(2));
t136 = -t124 * t102 + t126 * t125;
t90 = 0.1e1 / t136 ^ 2;
t158 = t88 * t90;
t122 = cos(pkin(6));
t117 = 0.1e1 / t122 ^ 2;
t121 = sin(pkin(6));
t115 = t121 ^ 2;
t120 = t126 ^ 2;
t107 = t120 * t115 * t117 + 0.1e1;
t119 = t124 ^ 2;
t150 = 0.1e1 / t107 ^ 2 * t119;
t157 = t117 * t150;
t147 = t126 * t121;
t106 = atan2(t147, t122);
t100 = sin(t106);
t101 = cos(t106);
t87 = t100 * t147 + t101 * t122;
t84 = 0.1e1 / t87;
t89 = 0.1e1 / t136;
t116 = 0.1e1 / t122;
t85 = 0.1e1 / t87 ^ 2;
t154 = t90 * t94;
t135 = t102 * qJD(2);
t143 = qJD(2) * t126;
t145 = qJD(1) * t126;
t146 = qJD(1) * t124;
t76 = -t103 * t145 + t123 * t146 + t124 * t135 - t125 * t143;
t139 = t76 * t154;
t93 = -t126 * t102 - t124 * t125;
t99 = t103 * qJD(2);
t77 = t93 * qJD(1) - t123 * t143 - t124 * t99;
t91 = t89 * t90;
t153 = t91 * t77;
t80 = 0.1e1 + t158;
t156 = (-t88 * t153 - t139) / t80 ^ 2;
t152 = t124 * t85;
t151 = t126 * t85;
t149 = t115 * t116;
t144 = qJD(2) * t124;
t142 = 0.2e1 * t93 * t94;
t104 = 0.1e1 / t107;
t138 = (t104 - 0.1e1) * t121;
t137 = -0.2e1 * t116 * t157;
t75 = (-t101 * t104 * t126 * t149 + t100 * t138) * t124;
t114 = t121 * t115;
t92 = t126 * t103 - t124 * t123;
t86 = t84 * t85;
t83 = t119 * t115 * t85 + 0.1e1;
t78 = 0.1e1 / t80;
t74 = qJD(1) * t75;
t1 = [(-t104 * t116 * t121 + t114 * t137) * t145, 0, 0, 0, 0, 0; (0.2e1 * (-t126 * t84 + t75 * t152) / t83 ^ 2 * (-t119 * t74 * t86 + t145 * t152) * t115 + ((0.2e1 * t124 * t75 * t86 - t151) * t74 + (-t75 * t151 + (-t84 + (-t114 * t157 - t138) * t100 * t151 - (t120 * t115 ^ 2 * t137 + (-t150 + (0.2e1 * t119 - t120) * t104) * t149) * t85 * t101) * t124) * qJD(1)) / t83) * t121, 0, 0, 0, 0, 0; (t90 * t142 - 0.2e1 * t92 * t89) * t156 + ((t93 * t76 - t92 * t77) * t90 + (-t103 * t146 - t123 * t145 - t125 * t144 - t126 * t135) * t89 - (-t136 * qJD(1) + t123 * t144 - t126 * t99) * t154 + t142 * t153) * t78, 0.2e1 * (-t136 * t89 - t158) * t156 + (-0.2e1 * t139 + (-t136 * t90 - 0.2e1 * t88 * t91 + t89) * t77) * t78, 0, 0, 0, 0;];
JaD_rot  = t1;
