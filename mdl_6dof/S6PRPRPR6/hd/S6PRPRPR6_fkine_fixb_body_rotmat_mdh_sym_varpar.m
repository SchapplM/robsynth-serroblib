% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRPRPR6 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:00
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PRPRPR6_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR6_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:00:41
	% EndTime: 2020-11-04 21:00:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:00:41
	% EndTime: 2020-11-04 21:00:41
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t79 = cos(pkin(10));
	t78 = sin(pkin(10));
	t1 = [t79, -t78, 0, 0; t78, t79, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:00:41
	% EndTime: 2020-11-04 21:00:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t80 = sin(pkin(10));
	t81 = sin(pkin(6));
	t89 = t80 * t81;
	t82 = cos(pkin(10));
	t88 = t82 * t81;
	t83 = cos(pkin(6));
	t84 = sin(qJ(2));
	t87 = t83 * t84;
	t85 = cos(qJ(2));
	t86 = t83 * t85;
	t1 = [-t80 * t87 + t82 * t85, -t80 * t86 - t82 * t84, t89, t82 * pkin(1) + pkin(7) * t89 + 0; t80 * t85 + t82 * t87, -t80 * t84 + t82 * t86, -t88, t80 * pkin(1) - pkin(7) * t88 + 0; t81 * t84, t81 * t85, t83, t83 * pkin(7) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:00:41
	% EndTime: 2020-11-04 21:00:41
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (23->23), mult. (43->34), div. (0->0), fcn. (56->6), ass. (0->15)
	t90 = sin(pkin(10));
	t103 = t90 * pkin(2);
	t92 = cos(pkin(10));
	t102 = t92 * pkin(2);
	t91 = sin(pkin(6));
	t101 = t90 * t91;
	t100 = t92 * t91;
	t93 = cos(pkin(6));
	t94 = sin(qJ(2));
	t99 = t93 * t94;
	t95 = cos(qJ(2));
	t98 = t93 * t95;
	t97 = t90 * qJ(3);
	t96 = t92 * qJ(3);
	t1 = [t101, t90 * t99 - t92 * t95, t90 * t98 + t92 * t94, (t93 * t97 + t102) * t95 + (-t93 * t103 + t96) * t94 + pkin(7) * t101 + t92 * pkin(1) + 0; -t100, -t90 * t95 - t92 * t99, t90 * t94 - t92 * t98, (-t93 * t96 + t103) * t95 + (t93 * t102 + t97) * t94 - pkin(7) * t100 + t90 * pkin(1) + 0; t93, -t91 * t94, -t91 * t95, t93 * pkin(7) + qJ(1) + 0 + (pkin(2) * t94 - qJ(3) * t95) * t91; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:00:41
	% EndTime: 2020-11-04 21:00:41
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (36->28), mult. (64->48), div. (0->0), fcn. (85->8), ass. (0->22)
	t106 = sin(pkin(10));
	t124 = t106 * qJ(3);
	t107 = sin(pkin(6));
	t110 = sin(qJ(4));
	t123 = t107 * t110;
	t112 = cos(qJ(4));
	t122 = t107 * t112;
	t113 = cos(qJ(2));
	t121 = t107 * t113;
	t114 = pkin(3) + pkin(7);
	t120 = t107 * t114;
	t108 = cos(pkin(10));
	t119 = t108 * qJ(3);
	t109 = cos(pkin(6));
	t111 = sin(qJ(2));
	t118 = t109 * t111;
	t117 = t109 * t113;
	t115 = pkin(2) + pkin(8);
	t116 = t109 * t115;
	t105 = t106 * t117 + t108 * t111;
	t104 = t106 * t111 - t108 * t117;
	t1 = [t105 * t110 + t106 * t122, t105 * t112 - t106 * t123, -t106 * t118 + t108 * t113, (t108 * t115 + t109 * t124) * t113 + (-t106 * t116 + t119) * t111 + t106 * t120 + t108 * pkin(1) + 0; t104 * t110 - t108 * t122, t104 * t112 + t108 * t123, t106 * t113 + t108 * t118, (t106 * t115 - t109 * t119) * t113 + (t108 * t116 + t124) * t111 - t108 * t120 + t106 * pkin(1) + 0; t109 * t112 - t110 * t121, -t109 * t110 - t112 * t121, t107 * t111, t114 * t109 + qJ(1) + 0 + (-qJ(3) * t113 + t111 * t115) * t107; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:00:41
	% EndTime: 2020-11-04 21:00:41
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (65->41), mult. (134->66), div. (0->0), fcn. (178->10), ass. (0->33)
	t136 = sin(pkin(10));
	t156 = t136 * qJ(3);
	t137 = sin(pkin(6));
	t141 = sin(qJ(4));
	t155 = t137 * t141;
	t142 = sin(qJ(2));
	t154 = t137 * t142;
	t143 = cos(qJ(4));
	t153 = t137 * t143;
	t144 = cos(qJ(2));
	t152 = t137 * t144;
	t145 = pkin(3) + pkin(7);
	t151 = t137 * t145;
	t139 = cos(pkin(10));
	t150 = t139 * qJ(3);
	t140 = cos(pkin(6));
	t149 = t140 * t142;
	t148 = t140 * t144;
	t146 = pkin(2) + pkin(8);
	t147 = t140 * t146;
	t138 = cos(pkin(11));
	t135 = sin(pkin(11));
	t134 = t140 * t143 - t141 * t152;
	t133 = t140 * t141 + t143 * t152;
	t132 = -t136 * t149 + t139 * t144;
	t131 = t136 * t148 + t139 * t142;
	t130 = t136 * t144 + t139 * t149;
	t129 = t136 * t142 - t139 * t148;
	t128 = t129 * t141 - t139 * t153;
	t127 = t129 * t143 + t139 * t155;
	t126 = t131 * t141 + t136 * t153;
	t125 = -t131 * t143 + t136 * t155;
	t1 = [t126 * t138 + t132 * t135, -t126 * t135 + t132 * t138, t125, t126 * pkin(4) + t125 * qJ(5) + (t139 * t146 + t140 * t156) * t144 + (-t136 * t147 + t150) * t142 + t136 * t151 + t139 * pkin(1) + 0; t128 * t138 + t130 * t135, -t128 * t135 + t130 * t138, -t127, t128 * pkin(4) - t127 * qJ(5) + (t136 * t146 - t140 * t150) * t144 + (t139 * t147 + t156) * t142 - t139 * t151 + t136 * pkin(1) + 0; t134 * t138 + t135 * t154, -t134 * t135 + t138 * t154, t133, t134 * pkin(4) + t133 * qJ(5) + t145 * t140 + qJ(1) + 0 + (-qJ(3) * t144 + t142 * t146) * t137; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:00:41
	% EndTime: 2020-11-04 21:00:41
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (88->47), mult. (150->70), div. (0->0), fcn. (196->12), ass. (0->37)
	t193 = sin(pkin(11)) * pkin(5);
	t172 = sin(pkin(10));
	t192 = t172 * qJ(3);
	t173 = sin(pkin(6));
	t177 = sin(qJ(4));
	t191 = t173 * t177;
	t178 = sin(qJ(2));
	t190 = t173 * t178;
	t179 = cos(qJ(4));
	t189 = t173 * t179;
	t180 = cos(qJ(2));
	t188 = t173 * t180;
	t181 = pkin(3) + pkin(7);
	t187 = t173 * t181;
	t174 = cos(pkin(10));
	t186 = t174 * qJ(3);
	t175 = cos(pkin(6));
	t185 = t175 * t178;
	t184 = t175 * t180;
	t182 = pkin(2) + pkin(8);
	t183 = t175 * t182;
	t176 = -pkin(9) - qJ(5);
	t170 = pkin(11) + qJ(6);
	t169 = cos(t170);
	t168 = sin(t170);
	t167 = cos(pkin(11)) * pkin(5) + pkin(4);
	t166 = t175 * t179 - t177 * t188;
	t165 = t175 * t177 + t179 * t188;
	t164 = -t172 * t185 + t174 * t180;
	t163 = t172 * t184 + t174 * t178;
	t162 = t172 * t180 + t174 * t185;
	t161 = t172 * t178 - t174 * t184;
	t160 = t161 * t177 - t174 * t189;
	t159 = t161 * t179 + t174 * t191;
	t158 = t163 * t177 + t172 * t189;
	t157 = -t163 * t179 + t172 * t191;
	t1 = [t158 * t169 + t164 * t168, -t158 * t168 + t164 * t169, t157, t158 * t167 - t157 * t176 + t164 * t193 + (t174 * t182 + t175 * t192) * t180 + (-t172 * t183 + t186) * t178 + t172 * t187 + t174 * pkin(1) + 0; t160 * t169 + t162 * t168, -t160 * t168 + t162 * t169, -t159, t160 * t167 + t159 * t176 + t162 * t193 + (t172 * t182 - t175 * t186) * t180 + (t174 * t183 + t192) * t178 - t174 * t187 + t172 * pkin(1) + 0; t166 * t169 + t168 * t190, -t166 * t168 + t169 * t190, t165, -t165 * t176 + t166 * t167 + t181 * t175 + qJ(1) + 0 + (-qJ(3) * t180 + (t182 + t193) * t178) * t173; 0, 0, 0, 1;];
	Tc_mdh = t1;
end