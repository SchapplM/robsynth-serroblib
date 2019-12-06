% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PPRRP2
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
%   Wie in S5PPRRP2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:09
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5PPRRP2_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP2_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP2_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRRP2_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP2_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:09:32
	% EndTime: 2019-12-05 15:09:32
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:09:32
	% EndTime: 2019-12-05 15:09:32
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:09:32
	% EndTime: 2019-12-05 15:09:32
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:09:32
	% EndTime: 2019-12-05 15:09:32
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:09:32
	% EndTime: 2019-12-05 15:09:32
	% DurationCPUTime: 0.34s
	% Computational Cost: add. (973->43), mult. (1166->102), div. (229->12), fcn. (1315->9), ass. (0->60)
	t104 = sin(pkin(7));
	t103 = pkin(8) + qJ(3);
	t100 = cos(t103);
	t99 = sin(t103);
	t95 = t99 ^ 2;
	t130 = t95 / t100 ^ 2;
	t120 = 0.1e1 + t130;
	t143 = t104 * t120;
	t142 = qJD(3) * t99;
	t106 = sin(qJ(4));
	t107 = cos(qJ(4));
	t105 = cos(pkin(7));
	t125 = t100 * t105;
	t84 = t104 * t106 + t107 * t125;
	t128 = t104 * t99;
	t87 = atan2(-t128, -t100);
	t85 = sin(t87);
	t86 = cos(t87);
	t70 = -t86 * t100 - t85 * t128;
	t67 = 0.1e1 / t70;
	t80 = 0.1e1 / t84;
	t101 = t104 ^ 2;
	t90 = t101 * t130 + 0.1e1;
	t88 = 0.1e1 / t90;
	t72 = t88 * t143;
	t117 = t104 * t72 - 0.1e1;
	t127 = t104 - t72;
	t129 = t100 * t85;
	t131 = t86 * t99;
	t61 = -t117 * t131 - t127 * t129;
	t140 = 0.2e1 * t61;
	t68 = 0.1e1 / t70 ^ 2;
	t81 = 0.1e1 / t84 ^ 2;
	t102 = t105 ^ 2;
	t123 = qJD(3) * t100;
	t135 = t68 * t99;
	t121 = t86 * t128;
	t71 = qJD(3) * t72;
	t60 = (-t121 + t129) * t71 + (-t104 * t129 + t131) * qJD(3);
	t138 = t60 * t67 * t68;
	t66 = t102 * t95 * t68 + 0.1e1;
	t139 = (t123 * t135 - t95 * t138) * t102 / t66 ^ 2;
	t83 = -t104 * t107 + t106 * t125;
	t132 = t81 * t83;
	t119 = t105 * t142;
	t126 = qJD(4) * t83;
	t77 = -t107 * t119 - t126;
	t133 = t77 * t80 * t81;
	t79 = t83 ^ 2;
	t75 = t79 * t81 + 0.1e1;
	t76 = -t84 * qJD(4) + t106 * t119;
	t137 = (-t76 * t132 - t79 * t133) / t75 ^ 2;
	t64 = 0.1e1 / t66;
	t136 = t64 * t68;
	t134 = t76 * t81;
	t122 = -0.2e1 * t137;
	t116 = -t106 * t80 + t107 * t132;
	t73 = 0.1e1 / t75;
	t62 = 0.2e1 * (t88 - t120 / t90 ^ 2 * t101) / t100 * t142 * t143;
	t1 = [0, 0, t62, 0, 0; 0, 0, ((-0.2e1 * t67 * t139 + (-qJD(3) * t61 - t60) * t136) * t100 + (t68 * t139 * t140 + (t62 * t68 * t121 + t138 * t140 - qJD(3) * t67 - (t127 * qJD(3) + t117 * t71) * t85 * t135) * t64 - (t62 * t85 + (qJD(3) + (-0.2e1 * t104 + t72) * t71) * t86) * t100 * t136) * t99) * t105, 0, 0; 0, 0, (t116 * t99 * t122 + (t116 * t123 + ((t77 - t126) * t81 * t106 + (-qJD(4) * t80 - 0.2e1 * t83 * t133 - t134) * t107) * t99) * t73) * t105, t122 + 0.2e1 * (-t73 * t134 + (-t73 * t133 - t81 * t137) * t83) * t83, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:09:32
	% EndTime: 2019-12-05 15:09:33
	% DurationCPUTime: 0.67s
	% Computational Cost: add. (2496->86), mult. (3936->202), div. (821->15), fcn. (5024->9), ass. (0->91)
	t134 = pkin(8) + qJ(3);
	t132 = cos(t134);
	t139 = cos(pkin(7));
	t141 = cos(qJ(4));
	t167 = t139 * t141;
	t138 = sin(pkin(7));
	t140 = sin(qJ(4));
	t170 = t138 * t140;
	t121 = t132 * t167 + t170;
	t115 = 0.1e1 / t121 ^ 2;
	t131 = sin(t134);
	t127 = t131 ^ 2;
	t133 = t139 ^ 2;
	t175 = t127 * t133;
	t110 = t115 * t175 + 0.1e1;
	t166 = qJD(3) * t132;
	t168 = t139 * t140;
	t169 = t138 * t141;
	t120 = t132 * t168 - t169;
	t164 = qJD(3) * t141;
	t154 = t131 * t164;
	t102 = -qJD(4) * t120 - t139 * t154;
	t114 = 0.1e1 / t121;
	t180 = t102 * t114 * t115;
	t190 = 0.1e1 / t110 ^ 2 * (t115 * t131 * t166 - t127 * t180) * t133;
	t128 = 0.1e1 / t131;
	t117 = t132 * t170 + t167;
	t119 = t132 * t169 - t168;
	t135 = 0.1e1 / t140;
	t136 = 0.1e1 / t140 ^ 2;
	t171 = t136 * t141;
	t150 = t117 * t171 - t119 * t135;
	t189 = t128 * t150;
	t112 = t117 ^ 2;
	t129 = 0.1e1 / t131 ^ 2;
	t174 = t129 * t136;
	t111 = t112 * t174 + 0.1e1;
	t130 = t128 / t127;
	t137 = t135 * t136;
	t158 = t117 * t174;
	t162 = qJD(4) * t141;
	t152 = t132 * t162;
	t173 = t131 * t138;
	t99 = -t138 * t152 + (qJD(3) * t173 + qJD(4) * t139) * t140;
	t188 = -0.2e1 / t111 ^ 2 * (-t99 * t158 + (-t129 * t137 * t162 - t130 * t136 * t166) * t112);
	t172 = t131 * t140;
	t109 = atan2(-t117, t172);
	t104 = cos(t109);
	t103 = sin(t109);
	t179 = t103 * t117;
	t98 = t172 * t104 - t179;
	t95 = 0.1e1 / t98;
	t187 = t139 * t131 * t95;
	t96 = 0.1e1 / t98 ^ 2;
	t186 = 0.2e1 * t120;
	t153 = t140 * t166;
	t149 = t131 * t162 + t153;
	t107 = 0.1e1 / t111;
	t87 = (t99 * t128 * t135 + t149 * t158) * t107;
	t183 = t117 * t87;
	t84 = (-t172 * t87 + t99) * t103 + (t149 - t183) * t104;
	t97 = t95 * t96;
	t185 = t84 * t97;
	t165 = qJD(3) * t139;
	t155 = t131 * t165;
	t163 = qJD(4) * t140;
	t101 = -t138 * t163 - t139 * t152 + t140 * t155;
	t184 = t101 * t96;
	t182 = t120 * t96;
	t156 = t129 * t132 * t135;
	t151 = t117 * t156 + t138;
	t94 = t151 * t107;
	t181 = t138 - t94;
	t178 = t103 * t131;
	t177 = t104 * t117;
	t176 = t104 * t132;
	t113 = t120 ^ 2;
	t92 = t113 * t96 + 0.1e1;
	t161 = 0.2e1 * (-t101 * t182 - t113 * t185) / t92 ^ 2;
	t160 = t131 * t186;
	t159 = t103 * t182;
	t157 = t141 * t175;
	t105 = 0.1e1 / t110;
	t100 = t117 * qJD(4) + t138 * t154;
	t90 = 0.1e1 / t92;
	t89 = t107 * t189;
	t86 = -t94 * t177 + (t178 * t181 + t176) * t140;
	t85 = (-t117 * t89 + t131 * t141) * t104 + (-t89 * t172 - t119) * t103;
	t83 = t151 * t188 + (-t99 * t156 + (-t152 * t174 + (-0.2e1 * t130 * t132 ^ 2 - t128) * t135 * qJD(3)) * t117) * t107;
	t81 = t188 * t189 + (-t150 * t129 * t166 + (-t99 * t171 + t100 * t135 + (t119 * t171 + (-0.2e1 * t137 * t141 ^ 2 - t135) * t117) * qJD(4)) * t128) * t107;
	t1 = [0, 0, t83, t81, 0; 0, 0, t86 * t161 * t182 + (t86 * t184 + (-(-t177 * t83 + (t104 * t99 + t179 * t87) * t94) * t96 + 0.2e1 * t86 * t185) * t120 + (-t187 - (t103 * t173 - t178 * t94 + t176) * t182) * t162) * t90 + (t161 * t187 + ((-t95 * t165 - (qJD(3) * t181 - t87) * t159) * t132 + (-(-t103 * t83 + (t181 * t87 - qJD(3)) * t104) * t182 + t139 * t96 * t84) * t131) * t90) * t140, (-t121 * t95 + t182 * t85) * t161 + (t85 * t184 + t102 * t95 + (t186 * t85 * t97 - t121 * t96) * t84 - (t132 * t164 - t117 * t81 - t119 * t87 + t89 * t99 + (-t87 * t89 - qJD(4)) * t172) * t104 * t182 - (t100 + (-t153 + t183) * t89 + (-t140 * t81 + (-qJD(4) * t89 - t87) * t141) * t131) * t159) * t90, 0; 0, 0, 0.2e1 * (t114 * t132 * t139 + t115 * t157) * t190 + (0.2e1 * t157 * t180 + t114 * t155 + (t163 * t175 + (t102 * t139 - 0.2e1 * t133 * t154) * t132) * t115) * t105, (t115 * t160 * t190 + (t160 * t180 + (t101 * t131 - t120 * t166) * t115) * t105) * t139, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end