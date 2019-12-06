% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PRPRR3
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S5PRPRR3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5PRPRR3_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR3_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR3_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPRR3_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR3_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:48:12
	% EndTime: 2019-12-05 15:48:12
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:48:12
	% EndTime: 2019-12-05 15:48:12
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:48:12
	% EndTime: 2019-12-05 15:48:12
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:48:12
	% EndTime: 2019-12-05 15:48:12
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:48:13
	% EndTime: 2019-12-05 15:48:13
	% DurationCPUTime: 0.35s
	% Computational Cost: add. (973->43), mult. (1166->102), div. (229->12), fcn. (1315->9), ass. (0->60)
	t105 = sin(pkin(8));
	t104 = qJ(2) + pkin(9);
	t101 = cos(t104);
	t100 = sin(t104);
	t96 = t100 ^ 2;
	t133 = t96 / t101 ^ 2;
	t121 = 0.1e1 + t133;
	t144 = t105 * t121;
	t143 = qJD(2) * t100;
	t107 = sin(qJ(4));
	t108 = cos(qJ(4));
	t106 = cos(pkin(8));
	t127 = t101 * t106;
	t85 = t105 * t107 + t108 * t127;
	t126 = t105 * t100;
	t88 = atan2(-t126, -t101);
	t86 = sin(t88);
	t87 = cos(t88);
	t71 = -t87 * t101 - t86 * t126;
	t68 = 0.1e1 / t71;
	t81 = 0.1e1 / t85;
	t102 = t105 ^ 2;
	t91 = t102 * t133 + 0.1e1;
	t89 = 0.1e1 / t91;
	t73 = t89 * t144;
	t119 = t105 * t73 - 0.1e1;
	t129 = t105 - t73;
	t130 = t87 * t100;
	t131 = t101 * t86;
	t62 = -t119 * t130 - t129 * t131;
	t141 = 0.2e1 * t62;
	t69 = 0.1e1 / t71 ^ 2;
	t82 = 0.1e1 / t85 ^ 2;
	t103 = t106 ^ 2;
	t124 = qJD(2) * t101;
	t132 = t100 * t69;
	t122 = t87 * t126;
	t72 = qJD(2) * t73;
	t61 = (-t122 + t131) * t72 + (-t105 * t131 + t130) * qJD(2);
	t139 = t61 * t68 * t69;
	t67 = t103 * t96 * t69 + 0.1e1;
	t140 = (t124 * t132 - t96 * t139) * t103 / t67 ^ 2;
	t84 = -t105 * t108 + t107 * t127;
	t134 = t82 * t84;
	t118 = t106 * t143;
	t128 = qJD(4) * t84;
	t78 = -t108 * t118 - t128;
	t135 = t78 * t81 * t82;
	t80 = t84 ^ 2;
	t76 = t80 * t82 + 0.1e1;
	t77 = -t85 * qJD(4) + t107 * t118;
	t138 = (-t77 * t134 - t80 * t135) / t76 ^ 2;
	t65 = 0.1e1 / t67;
	t137 = t65 * t69;
	t136 = t77 * t82;
	t123 = -0.2e1 * t138;
	t117 = -t107 * t81 + t108 * t134;
	t74 = 0.1e1 / t76;
	t63 = 0.2e1 * (t89 - t121 / t91 ^ 2 * t102) / t101 * t143 * t144;
	t1 = [0, t63, 0, 0, 0; 0, ((-0.2e1 * t68 * t140 + (-qJD(2) * t62 - t61) * t137) * t101 + (t69 * t140 * t141 + (t63 * t69 * t122 + t139 * t141 - qJD(2) * t68 - (t129 * qJD(2) + t119 * t72) * t86 * t132) * t65 - (t63 * t86 + (qJD(2) + (-0.2e1 * t105 + t73) * t72) * t87) * t101 * t137) * t100) * t106, 0, 0, 0; 0, (t117 * t74 * t124 + (t117 * t123 + ((t78 - t128) * t82 * t107 + (-qJD(4) * t81 - 0.2e1 * t84 * t135 - t136) * t108) * t74) * t100) * t106, 0, t123 + 0.2e1 * (-t74 * t136 + (-t74 * t135 - t82 * t138) * t84) * t84, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:48:13
	% EndTime: 2019-12-05 15:48:13
	% DurationCPUTime: 0.33s
	% Computational Cost: add. (1374->45), mult. (1346->104), div. (247->12), fcn. (1511->9), ass. (0->63)
	t138 = sin(pkin(8));
	t136 = qJ(2) + pkin(9);
	t129 = sin(t136);
	t125 = t129 ^ 2;
	t130 = cos(t136);
	t163 = t125 / t130 ^ 2;
	t151 = 0.1e1 + t163;
	t175 = t138 * t151;
	t174 = qJD(2) * t129;
	t137 = qJ(4) + qJ(5);
	t131 = sin(t137);
	t132 = cos(t137);
	t139 = cos(pkin(8));
	t162 = t130 * t139;
	t114 = t138 * t131 + t132 * t162;
	t133 = t138 ^ 2;
	t121 = t133 * t163 + 0.1e1;
	t119 = 0.1e1 / t121;
	t102 = t119 * t175;
	t161 = t138 * t129;
	t117 = atan2(-t161, -t130);
	t116 = cos(t117);
	t115 = sin(t117);
	t164 = t115 * t130;
	t148 = t116 * t129 - t138 * t164;
	t156 = t116 * t161;
	t149 = -t156 + t164;
	t91 = t149 * t102 + t148;
	t172 = 0.2e1 * t91;
	t100 = -t115 * t161 - t116 * t130;
	t97 = 0.1e1 / t100;
	t110 = 0.1e1 / t114;
	t98 = 0.1e1 / t100 ^ 2;
	t111 = 0.1e1 / t114 ^ 2;
	t134 = t139 ^ 2;
	t158 = qJD(2) * t130;
	t167 = t129 * t98;
	t101 = qJD(2) * t102;
	t90 = t148 * qJD(2) + t149 * t101;
	t170 = t90 * t97 * t98;
	t96 = t134 * t125 * t98 + 0.1e1;
	t171 = (-t125 * t170 + t158 * t167) * t134 / t96 ^ 2;
	t94 = 0.1e1 / t96;
	t169 = t94 * t98;
	t155 = t131 * t162;
	t113 = -t138 * t132 + t155;
	t109 = t113 ^ 2;
	t105 = t109 * t111 + 0.1e1;
	t135 = qJD(4) + qJD(5);
	t153 = t139 * t174;
	t106 = -t114 * t135 + t131 * t153;
	t165 = t111 * t113;
	t107 = -t135 * t155 + (t135 * t138 - t153) * t132;
	t166 = t107 * t110 * t111;
	t168 = 0.1e1 / t105 ^ 2 * (-t106 * t165 - t109 * t166);
	t159 = t102 - t138;
	t157 = -0.2e1 * t168;
	t152 = t102 * t138 - 0.1e1;
	t150 = -t110 * t131 + t132 * t165;
	t103 = 0.1e1 / t105;
	t93 = 0.2e1 * (t119 - t151 / t121 ^ 2 * t133) / t130 * t174 * t175;
	t88 = t157 + 0.2e1 * (-t103 * t106 * t111 + (-t103 * t166 - t111 * t168) * t113) * t113;
	t1 = [0, t93, 0, 0, 0; 0, ((-0.2e1 * t97 * t171 + (-qJD(2) * t91 - t90) * t169) * t130 + (t98 * t171 * t172 + (t93 * t98 * t156 + t170 * t172 - qJD(2) * t97 - (-t159 * qJD(2) + t152 * t101) * t115 * t167) * t94 - (t115 * t93 + (-t152 * qJD(2) + t159 * t101) * t116) * t130 * t169) * t129) * t139, 0, 0, 0; 0, (t150 * t129 * t157 + (t150 * t158 + ((-t110 * t135 - 0.2e1 * t113 * t166) * t132 + (-t106 * t132 + (-t113 * t135 + t107) * t131) * t111) * t129) * t103) * t139, 0, t88, t88;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end