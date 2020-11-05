% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRPRPR4 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:59
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PRPRPR4_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR4_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR4_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR4_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:59:56
	% EndTime: 2020-11-04 20:59:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:59:56
	% EndTime: 2020-11-04 20:59:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t83 = cos(pkin(10));
	t82 = sin(pkin(10));
	t1 = [t83, -t82, 0, 0; t82, t83, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:59:56
	% EndTime: 2020-11-04 20:59:56
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t84 = sin(pkin(10));
	t85 = sin(pkin(6));
	t93 = t84 * t85;
	t86 = cos(pkin(10));
	t92 = t86 * t85;
	t87 = cos(pkin(6));
	t88 = sin(qJ(2));
	t91 = t87 * t88;
	t89 = cos(qJ(2));
	t90 = t87 * t89;
	t1 = [-t84 * t91 + t86 * t89, -t84 * t90 - t86 * t88, t93, t86 * pkin(1) + pkin(7) * t93 + 0; t84 * t89 + t86 * t91, -t84 * t88 + t86 * t90, -t92, t84 * pkin(1) - pkin(7) * t92 + 0; t85 * t88, t85 * t89, t87, t87 * pkin(7) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:59:56
	% EndTime: 2020-11-04 20:59:56
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (29->27), mult. (64->47), div. (0->0), fcn. (85->8), ass. (0->19)
	t97 = sin(pkin(10));
	t111 = t97 * pkin(2);
	t98 = sin(pkin(6));
	t110 = t97 * t98;
	t100 = cos(pkin(10));
	t109 = t100 * t98;
	t102 = sin(qJ(2));
	t108 = t102 * t98;
	t107 = t97 * qJ(3);
	t101 = cos(pkin(6));
	t106 = t100 * t101;
	t105 = t101 * t102;
	t103 = cos(qJ(2));
	t104 = t101 * t103;
	t99 = cos(pkin(11));
	t96 = sin(pkin(11));
	t95 = t100 * t105 + t97 * t103;
	t94 = -t100 * t103 + t97 * t105;
	t1 = [t96 * t110 - t94 * t99, t99 * t110 + t94 * t96, t100 * t102 + t97 * t104, (t100 * pkin(2) + t101 * t107) * t103 + (t100 * qJ(3) - t101 * t111) * t102 + pkin(7) * t110 + t100 * pkin(1) + 0; -t96 * t109 + t95 * t99, -t99 * t109 - t95 * t96, -t100 * t104 + t97 * t102, (-qJ(3) * t106 + t111) * t103 + (pkin(2) * t106 + t107) * t102 - pkin(7) * t109 + t97 * pkin(1) + 0; t101 * t96 + t99 * t108, t101 * t99 - t96 * t108, -t98 * t103, t101 * pkin(7) + qJ(1) + 0 + (pkin(2) * t102 - qJ(3) * t103) * t98; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:59:56
	% EndTime: 2020-11-04 20:59:56
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (54->31), mult. (72->49), div. (0->0), fcn. (93->10), ass. (0->23)
	t115 = cos(pkin(11)) * pkin(3) + pkin(2);
	t119 = sin(pkin(10));
	t133 = t119 * t115;
	t120 = sin(pkin(6));
	t132 = t119 * t120;
	t121 = cos(pkin(10));
	t131 = t120 * t121;
	t124 = sin(qJ(2));
	t130 = t120 * t124;
	t129 = t121 * t115;
	t122 = cos(pkin(6));
	t123 = qJ(3) + pkin(8);
	t128 = t122 * t123;
	t127 = t122 * t124;
	t125 = cos(qJ(2));
	t126 = t122 * t125;
	t118 = pkin(11) + qJ(4);
	t117 = cos(t118);
	t116 = sin(t118);
	t114 = sin(pkin(11)) * pkin(3) + pkin(7);
	t113 = t119 * t125 + t121 * t127;
	t112 = t119 * t127 - t121 * t125;
	t1 = [-t112 * t117 + t116 * t132, t112 * t116 + t117 * t132, t119 * t126 + t121 * t124, (t119 * t128 + t129) * t125 + (t121 * t123 - t122 * t133) * t124 + t114 * t132 + t121 * pkin(1) + 0; t113 * t117 - t116 * t131, -t113 * t116 - t117 * t131, t119 * t124 - t121 * t126, (-t121 * t128 + t133) * t125 + (t119 * t123 + t122 * t129) * t124 - t114 * t131 + t119 * pkin(1) + 0; t122 * t116 + t117 * t130, -t116 * t130 + t122 * t117, -t120 * t125, t114 * t122 + qJ(1) + 0 + (t115 * t124 - t123 * t125) * t120; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:59:56
	% EndTime: 2020-11-04 20:59:56
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (100->43), mult. (142->67), div. (0->0), fcn. (186->12), ass. (0->34)
	t145 = cos(pkin(11)) * pkin(3) + pkin(2);
	t150 = sin(pkin(10));
	t166 = t150 * t145;
	t151 = sin(pkin(6));
	t165 = t150 * t151;
	t153 = cos(pkin(10));
	t164 = t151 * t153;
	t156 = sin(qJ(2));
	t163 = t151 * t156;
	t157 = cos(qJ(2));
	t162 = t151 * t157;
	t161 = t153 * t145;
	t154 = cos(pkin(6));
	t155 = qJ(3) + pkin(8);
	t160 = t154 * t155;
	t159 = t154 * t156;
	t158 = t154 * t157;
	t152 = cos(pkin(12));
	t149 = sin(pkin(12));
	t148 = pkin(11) + qJ(4);
	t147 = cos(t148);
	t146 = sin(t148);
	t144 = sin(pkin(11)) * pkin(3) + pkin(7);
	t143 = t150 * t158 + t153 * t156;
	t142 = t150 * t157 + t153 * t159;
	t141 = t150 * t156 - t153 * t158;
	t140 = t150 * t159 - t153 * t157;
	t139 = t154 * t146 + t147 * t163;
	t138 = t146 * t163 - t154 * t147;
	t137 = -t140 * t147 + t146 * t165;
	t136 = t142 * t147 - t146 * t164;
	t135 = t142 * t146 + t147 * t164;
	t134 = t140 * t146 + t147 * t165;
	t1 = [t137 * t152 + t143 * t149, -t137 * t149 + t143 * t152, -t134, t137 * pkin(4) - t134 * qJ(5) + (t150 * t160 + t161) * t157 + (t153 * t155 - t154 * t166) * t156 + t144 * t165 + t153 * pkin(1) + 0; t136 * t152 + t141 * t149, -t136 * t149 + t141 * t152, t135, t136 * pkin(4) + t135 * qJ(5) + (-t153 * t160 + t166) * t157 + (t150 * t155 + t154 * t161) * t156 - t144 * t164 + t150 * pkin(1) + 0; t139 * t152 - t149 * t162, -t139 * t149 - t152 * t162, t138, t139 * pkin(4) + t138 * qJ(5) + t144 * t154 + qJ(1) + 0 + (t145 * t156 - t155 * t157) * t151; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:59:56
	% EndTime: 2020-11-04 20:59:57
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (123->49), mult. (158->71), div. (0->0), fcn. (204->14), ass. (0->38)
	t204 = sin(pkin(12)) * pkin(5);
	t179 = cos(pkin(11)) * pkin(3) + pkin(2);
	t187 = sin(pkin(10));
	t203 = t187 * t179;
	t188 = sin(pkin(6));
	t202 = t187 * t188;
	t189 = cos(pkin(10));
	t201 = t188 * t189;
	t193 = sin(qJ(2));
	t200 = t188 * t193;
	t194 = cos(qJ(2));
	t199 = t188 * t194;
	t198 = t189 * t179;
	t190 = cos(pkin(6));
	t191 = qJ(3) + pkin(8);
	t197 = t190 * t191;
	t196 = t190 * t193;
	t195 = t190 * t194;
	t192 = -pkin(9) - qJ(5);
	t185 = pkin(11) + qJ(4);
	t184 = pkin(12) + qJ(6);
	t183 = cos(t185);
	t182 = cos(t184);
	t181 = sin(t185);
	t180 = sin(t184);
	t178 = cos(pkin(12)) * pkin(5) + pkin(4);
	t177 = sin(pkin(11)) * pkin(3) + pkin(7);
	t176 = t187 * t195 + t189 * t193;
	t175 = t187 * t194 + t189 * t196;
	t174 = t187 * t193 - t189 * t195;
	t173 = t187 * t196 - t189 * t194;
	t172 = t190 * t181 + t183 * t200;
	t171 = t181 * t200 - t190 * t183;
	t170 = -t173 * t183 + t181 * t202;
	t169 = t175 * t183 - t181 * t201;
	t168 = t175 * t181 + t183 * t201;
	t167 = t173 * t181 + t183 * t202;
	t1 = [t170 * t182 + t176 * t180, -t170 * t180 + t176 * t182, -t167, t170 * t178 + t167 * t192 + t176 * t204 + (t187 * t197 + t198) * t194 + (t189 * t191 - t190 * t203) * t193 + t177 * t202 + t189 * pkin(1) + 0; t169 * t182 + t174 * t180, -t169 * t180 + t174 * t182, t168, t169 * t178 - t168 * t192 + t174 * t204 + (-t189 * t197 + t203) * t194 + (t187 * t191 + t190 * t198) * t193 - t177 * t201 + t187 * pkin(1) + 0; t172 * t182 - t180 * t199, -t172 * t180 - t182 * t199, t171, -t171 * t192 + t172 * t178 + t177 * t190 + qJ(1) + 0 + (t179 * t193 + (-t191 - t204) * t194) * t188; 0, 0, 0, 1;];
	Tc_mdh = t1;
end