% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PRPRP4
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
%   Wie in S5PRPRP4_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:24
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5PRPRP4_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP4_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP4_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPRP4_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP4_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:24:26
	% EndTime: 2019-10-24 10:24:26
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:24:26
	% EndTime: 2019-10-24 10:24:26
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:24:26
	% EndTime: 2019-10-24 10:24:26
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:24:25
	% EndTime: 2019-10-24 10:24:26
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:24:26
	% EndTime: 2019-10-24 10:24:26
	% DurationCPUTime: 0.36s
	% Computational Cost: add. (973->43), mult. (1166->102), div. (229->12), fcn. (1315->9), ass. (0->60)
	t105 = sin(pkin(7));
	t104 = qJ(2) + pkin(8);
	t101 = cos(t104);
	t100 = sin(t104);
	t96 = t100 ^ 2;
	t133 = t96 / t101 ^ 2;
	t121 = 0.1e1 + t133;
	t144 = t105 * t121;
	t143 = qJD(2) * t100;
	t107 = sin(qJ(4));
	t108 = cos(qJ(4));
	t106 = cos(pkin(7));
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
	% StartTime: 2019-10-24 10:24:26
	% EndTime: 2019-10-24 10:24:27
	% DurationCPUTime: 0.66s
	% Computational Cost: add. (2496->86), mult. (3936->202), div. (821->15), fcn. (5024->9), ass. (0->91)
	t135 = qJ(2) + pkin(8);
	t133 = cos(t135);
	t140 = cos(pkin(7));
	t142 = cos(qJ(4));
	t168 = t140 * t142;
	t139 = sin(pkin(7));
	t141 = sin(qJ(4));
	t171 = t139 * t141;
	t122 = t133 * t168 + t171;
	t116 = 0.1e1 / t122 ^ 2;
	t132 = sin(t135);
	t128 = t132 ^ 2;
	t134 = t140 ^ 2;
	t176 = t128 * t134;
	t111 = t116 * t176 + 0.1e1;
	t167 = qJD(2) * t133;
	t169 = t140 * t141;
	t170 = t139 * t142;
	t121 = t133 * t169 - t170;
	t165 = qJD(2) * t142;
	t155 = t132 * t165;
	t103 = -t121 * qJD(4) - t140 * t155;
	t115 = 0.1e1 / t122;
	t181 = t103 * t115 * t116;
	t191 = 0.1e1 / t111 ^ 2 * (t116 * t132 * t167 - t128 * t181) * t134;
	t129 = 0.1e1 / t132;
	t118 = t133 * t171 + t168;
	t120 = t133 * t170 - t169;
	t136 = 0.1e1 / t141;
	t137 = 0.1e1 / t141 ^ 2;
	t172 = t137 * t142;
	t151 = t118 * t172 - t120 * t136;
	t190 = t129 * t151;
	t163 = qJD(4) * t142;
	t153 = t133 * t163;
	t174 = t132 * t139;
	t100 = -t139 * t153 + (qJD(2) * t174 + qJD(4) * t140) * t141;
	t113 = t118 ^ 2;
	t130 = 0.1e1 / t132 ^ 2;
	t175 = t130 * t137;
	t112 = t113 * t175 + 0.1e1;
	t131 = t129 / t128;
	t138 = t136 * t137;
	t159 = t118 * t175;
	t189 = -0.2e1 / t112 ^ 2 * (-t100 * t159 + (-t130 * t138 * t163 - t131 * t137 * t167) * t113);
	t173 = t132 * t141;
	t110 = atan2(-t118, t173);
	t105 = cos(t110);
	t104 = sin(t110);
	t180 = t104 * t118;
	t99 = t105 * t173 - t180;
	t96 = 0.1e1 / t99;
	t188 = t140 * t132 * t96;
	t97 = 0.1e1 / t99 ^ 2;
	t187 = 0.2e1 * t121;
	t154 = t141 * t167;
	t150 = t132 * t163 + t154;
	t108 = 0.1e1 / t112;
	t88 = (t100 * t129 * t136 + t150 * t159) * t108;
	t184 = t118 * t88;
	t85 = (-t88 * t173 + t100) * t104 + (t150 - t184) * t105;
	t98 = t96 * t97;
	t186 = t85 * t98;
	t166 = qJD(2) * t140;
	t156 = t132 * t166;
	t164 = qJD(4) * t141;
	t102 = -t139 * t164 - t140 * t153 + t141 * t156;
	t185 = t102 * t97;
	t183 = t121 * t97;
	t157 = t130 * t133 * t136;
	t152 = t118 * t157 + t139;
	t95 = t152 * t108;
	t182 = t139 - t95;
	t179 = t104 * t132;
	t178 = t105 * t118;
	t177 = t105 * t133;
	t114 = t121 ^ 2;
	t93 = t114 * t97 + 0.1e1;
	t162 = 0.2e1 * (-t102 * t183 - t114 * t186) / t93 ^ 2;
	t161 = t132 * t187;
	t160 = t104 * t183;
	t158 = t142 * t176;
	t106 = 0.1e1 / t111;
	t101 = qJD(4) * t118 + t139 * t155;
	t91 = 0.1e1 / t93;
	t90 = t108 * t190;
	t87 = -t95 * t178 + (t182 * t179 + t177) * t141;
	t86 = (-t118 * t90 + t132 * t142) * t105 + (-t90 * t173 - t120) * t104;
	t84 = t152 * t189 + (-t100 * t157 + (-t153 * t175 + (-0.2e1 * t131 * t133 ^ 2 - t129) * t136 * qJD(2)) * t118) * t108;
	t82 = t189 * t190 + (-t151 * t130 * t167 + (-t100 * t172 + t101 * t136 + (t120 * t172 + (-0.2e1 * t138 * t142 ^ 2 - t136) * t118) * qJD(4)) * t129) * t108;
	t1 = [0, t84, 0, t82, 0; 0, t87 * t162 * t183 + (t87 * t185 + (-(-t84 * t178 + (t100 * t105 + t180 * t88) * t95) * t97 + 0.2e1 * t87 * t186) * t121 + (-t188 - (t104 * t174 - t95 * t179 + t177) * t183) * t163) * t91 + (t162 * t188 + ((-t96 * t166 - (t182 * qJD(2) - t88) * t160) * t133 + (-(-t104 * t84 + (t182 * t88 - qJD(2)) * t105) * t183 + t140 * t97 * t85) * t132) * t91) * t141, 0, (-t122 * t96 + t86 * t183) * t162 + (t86 * t185 + t103 * t96 + (t86 * t98 * t187 - t122 * t97) * t85 - (t133 * t165 + t100 * t90 - t118 * t82 - t120 * t88 + (-t88 * t90 - qJD(4)) * t173) * t105 * t183 - (t101 + (-t154 + t184) * t90 + (-t141 * t82 + (-qJD(4) * t90 - t88) * t142) * t132) * t160) * t91, 0; 0, 0.2e1 * (t115 * t133 * t140 + t116 * t158) * t191 + (0.2e1 * t158 * t181 + t115 * t156 + (t164 * t176 + (t103 * t140 - 0.2e1 * t134 * t155) * t133) * t116) * t106, 0, (t116 * t161 * t191 + (t161 * t181 + (t102 * t132 - t121 * t167) * t116) * t106) * t140, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end