% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR10
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 11:27
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRRR10_jacobia_transl_6_floatb_twist_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_jacobia_transl_6_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR10_jacobia_transl_6_floatb_twist_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_jacobia_transl_6_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 11:27:17
% EndTime: 2018-11-23 11:27:18
% DurationCPUTime: 1.50s
% Computational Cost: add. (6634->201), mult. (6681->307), div. (0->0), fcn. (6693->30), ass. (0->122)
t106 = sin(qJ(6));
t112 = cos(qJ(6));
t108 = sin(qJ(4));
t143 = pkin(8) + qJ(4);
t130 = sin(t143) / 0.2e1;
t144 = pkin(8) - qJ(4);
t136 = sin(t144);
t121 = t130 + t136 / 0.2e1;
t133 = cos(t143) / 0.2e1;
t139 = cos(t144);
t122 = t139 / 0.2e1 + t133;
t102 = sin(pkin(7));
t117 = cos(qJ(1));
t103 = sin(pkin(6));
t154 = cos(pkin(7));
t142 = t103 * t154;
t110 = sin(qJ(2));
t111 = sin(qJ(1));
t147 = pkin(6) + qJ(2);
t135 = cos(t147) / 0.2e1;
t148 = pkin(6) - qJ(2);
t141 = cos(t148);
t123 = t141 / 0.2e1 + t135;
t76 = t111 * t110 - t117 * t123;
t126 = t76 * t102 - t117 * t142;
t115 = cos(qJ(3));
t149 = t103 * t117;
t116 = cos(qJ(2));
t132 = sin(t147) / 0.2e1;
t138 = sin(t148);
t87 = t132 - t138 / 0.2e1;
t77 = t111 * t116 + t117 * t87;
t145 = pkin(7) + qJ(3);
t131 = sin(t145) / 0.2e1;
t146 = pkin(7) - qJ(3);
t137 = sin(t146);
t85 = t131 - t137 / 0.2e1;
t134 = cos(t145) / 0.2e1;
t140 = cos(t146);
t89 = t134 - t140 / 0.2e1;
t55 = t77 * t115 + t149 * t89 - t76 * t85;
t109 = sin(qJ(3));
t84 = t131 + t137 / 0.2e1;
t90 = t140 / 0.2e1 + t134;
t57 = t77 * t109 + t149 * t84 + t76 * t90;
t19 = t55 * t108 - t121 * t126 + t122 * t57;
t107 = sin(qJ(5));
t113 = cos(qJ(5));
t114 = cos(qJ(4));
t83 = t130 - t136 / 0.2e1;
t88 = t133 - t139 / 0.2e1;
t20 = t55 * t114 - t126 * t88 - t57 * t83;
t101 = sin(pkin(8));
t104 = cos(pkin(8));
t47 = t57 * t101 + t104 * t126;
t4 = t47 * t107 + t20 * t113;
t164 = t4 * t106 - t19 * t112;
t163 = -t19 * t106 - t4 * t112;
t160 = -t20 * t107 + t47 * t113;
t159 = r_i_i_C(3) + pkin(14);
t158 = pkin(12) * t101;
t157 = t102 * pkin(11);
t156 = t102 * t88;
t91 = t135 - t141 / 0.2e1;
t155 = t102 * t91;
t153 = t101 * t107;
t152 = t101 * t113;
t151 = t102 * t104;
t150 = t103 * t111;
t81 = t111 * t87 - t117 * t116;
t129 = t112 * r_i_i_C(1) - t106 * r_i_i_C(2) + pkin(5);
t128 = t106 * r_i_i_C(1) + t112 * r_i_i_C(2) + pkin(13);
t63 = t76 * t109 - t77 * t90;
t49 = -t63 * t101 + t151 * t77;
t79 = -t117 * t110 - t111 * t123;
t65 = -t79 * t109 + t81 * t90;
t50 = -t65 * t101 - t151 * t81;
t86 = t132 + t138 / 0.2e1;
t72 = -t86 * t109 + t91 * t90;
t62 = -t72 * t101 - t151 * t91;
t105 = cos(pkin(6));
t127 = t86 * t102 - t105 * t154;
t69 = t105 * t89 + t91 * t115 - t86 * t85;
t125 = t79 * t102 - t111 * t142;
t124 = t109 * t81 + t150 * t84 + t79 * t90;
t60 = t115 * t81 + t150 * t89 - t79 * t85;
t120 = t102 * t121;
t119 = -t159 * t107 - t129 * t113 - pkin(4);
t67 = t105 * t84 + t91 * t109 + t86 * t90;
t40 = -t114 * t69 + t127 * t88 + t67 * t83;
t118 = -t101 * t124 - t104 * t125;
t25 = -t114 * t60 + t124 * t83 + t125 * t88;
t73 = t86 * t115 + t91 * t85;
t66 = t79 * t115 + t81 * t85;
t64 = -t76 * t115 - t77 * t85;
t53 = -t67 * t101 - t104 * t127;
t46 = t73 * t114 + t155 * t88 + t72 * t83;
t45 = t73 * t108 + t120 * t91 - t122 * t72;
t43 = t67 * t114 + t69 * t83;
t42 = t67 * t108 - t122 * t69;
t39 = -t108 * t69 + t121 * t127 - t122 * t67;
t38 = t114 * t124 + t60 * t83;
t37 = t108 * t124 - t122 * t60;
t36 = -t114 * t57 - t55 * t83;
t35 = -t108 * t57 + t122 * t55;
t34 = t66 * t114 + t156 * t81 + t65 * t83;
t33 = t66 * t108 + t120 * t81 - t122 * t65;
t32 = t64 * t114 - t156 * t77 + t63 * t83;
t31 = t64 * t108 - t120 * t77 - t122 * t63;
t30 = t43 * t113 - t153 * t69;
t28 = t62 * t107 + t46 * t113;
t24 = -t108 * t60 + t121 * t125 - t122 * t124;
t18 = t53 * t107 + t40 * t113;
t16 = t38 * t113 - t153 * t60;
t14 = t36 * t113 + t153 * t55;
t12 = t50 * t107 + t34 * t113;
t10 = t49 * t107 + t32 * t113;
t8 = t107 * t118 + t25 * t113;
t7 = t25 * t107 - t113 * t118;
t2 = t24 * t106 + t8 * t112;
t1 = -t8 * t106 + t24 * t112;
t3 = [-t111 * pkin(1) - t77 * pkin(2) - t55 * pkin(3) - t20 * pkin(4) - t4 * pkin(5) + pkin(10) * t149 - t126 * pkin(11) - t47 * pkin(12) - t19 * pkin(13) + t163 * r_i_i_C(1) + t164 * r_i_i_C(2) + t159 * t160 (t33 * t106 + t12 * t112) * r_i_i_C(1) + (-t12 * t106 + t33 * t112) * r_i_i_C(2) + t12 * pkin(5) + t34 * pkin(4) + t33 * pkin(13) + t66 * pkin(3) + t79 * pkin(2) - t81 * t157 + t159 * (t34 * t107 - t50 * t113) + t50 * pkin(12) (t37 * t106 + t16 * t112) * r_i_i_C(1) + (-t16 * t106 + t37 * t112) * r_i_i_C(2) + t16 * pkin(5) + t38 * pkin(4) + t37 * pkin(13) + t124 * pkin(3) - t60 * t158 + t159 * (t38 * t107 + t152 * t60) t119 * t24 + t128 * t25, -t129 * t7 + t159 * t8, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t117 * pkin(1) - pkin(2) * t81 - pkin(3) * t60 + t25 * pkin(4) + t8 * pkin(5) + pkin(10) * t150 - t125 * pkin(11) + t118 * pkin(12) + t24 * pkin(13) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t159 * t7 (t10 * t112 + t31 * t106) * r_i_i_C(1) + (-t10 * t106 + t31 * t112) * r_i_i_C(2) + t10 * pkin(5) + t32 * pkin(4) + t31 * pkin(13) + t64 * pkin(3) - t76 * pkin(2) + t77 * t157 + t159 * (t32 * t107 - t49 * t113) + t49 * pkin(12) (t35 * t106 + t14 * t112) * r_i_i_C(1) + (-t14 * t106 + t35 * t112) * r_i_i_C(2) + t14 * pkin(5) + t36 * pkin(4) + t35 * pkin(13) - t57 * pkin(3) + t55 * t158 + t159 * (t36 * t107 - t152 * t55) t119 * t19 + t128 * t20, t129 * t160 + t159 * t4, -t164 * r_i_i_C(1) + t163 * r_i_i_C(2); 0 (t45 * t106 + t28 * t112) * r_i_i_C(1) + (-t28 * t106 + t45 * t112) * r_i_i_C(2) + t28 * pkin(5) + t46 * pkin(4) + t45 * pkin(13) + t73 * pkin(3) + t86 * pkin(2) - pkin(11) * t155 + t159 * (t46 * t107 - t62 * t113) + t62 * pkin(12) (t42 * t106 + t30 * t112) * r_i_i_C(1) + (-t30 * t106 + t42 * t112) * r_i_i_C(2) + t30 * pkin(5) + t43 * pkin(4) + t42 * pkin(13) + t67 * pkin(3) - t69 * t158 + t159 * (t43 * t107 + t152 * t69) t119 * t39 + t128 * t40, t159 * t18 + t129 * (-t40 * t107 + t53 * t113) (-t18 * t106 + t39 * t112) * r_i_i_C(1) + (-t39 * t106 - t18 * t112) * r_i_i_C(2);];
Ja_transl  = t3;
