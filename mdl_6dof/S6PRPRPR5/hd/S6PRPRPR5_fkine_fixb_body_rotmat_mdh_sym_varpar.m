% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRPRPR5 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:00
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PRPRPR5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:00:19
	% EndTime: 2020-11-04 21:00:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:00:19
	% EndTime: 2020-11-04 21:00:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t81 = cos(pkin(10));
	t80 = sin(pkin(10));
	t1 = [t81, -t80, 0, 0; t80, t81, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:00:19
	% EndTime: 2020-11-04 21:00:19
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t82 = sin(pkin(10));
	t83 = sin(pkin(6));
	t91 = t82 * t83;
	t84 = cos(pkin(10));
	t90 = t84 * t83;
	t85 = cos(pkin(6));
	t86 = sin(qJ(2));
	t89 = t85 * t86;
	t87 = cos(qJ(2));
	t88 = t85 * t87;
	t1 = [-t82 * t89 + t84 * t87, -t82 * t88 - t84 * t86, t91, t84 * pkin(1) + pkin(7) * t91 + 0; t82 * t87 + t84 * t89, -t82 * t86 + t84 * t88, -t90, t82 * pkin(1) - pkin(7) * t90 + 0; t83 * t86, t83 * t87, t85, t85 * pkin(7) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:00:19
	% EndTime: 2020-11-04 21:00:19
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (29->27), mult. (64->44), div. (0->0), fcn. (85->8), ass. (0->22)
	t95 = sin(pkin(10));
	t112 = t95 * pkin(2);
	t98 = cos(pkin(10));
	t111 = t98 * pkin(2);
	t96 = sin(pkin(6));
	t110 = t95 * t96;
	t109 = t96 * t98;
	t100 = sin(qJ(2));
	t108 = t100 * t96;
	t107 = t95 * qJ(3);
	t106 = t95 * t100;
	t101 = cos(qJ(2));
	t105 = t95 * t101;
	t104 = t98 * qJ(3);
	t103 = t98 * t100;
	t102 = t98 * t101;
	t99 = cos(pkin(6));
	t97 = cos(pkin(11));
	t94 = sin(pkin(11));
	t93 = -t99 * t106 + t102;
	t92 = t99 * t103 + t105;
	t1 = [t94 * t110 + t93 * t97, t97 * t110 - t93 * t94, t99 * t105 + t103, (t99 * t107 + t111) * t101 + (-t99 * t112 + t104) * t100 + pkin(7) * t110 + t98 * pkin(1) + 0; -t94 * t109 + t92 * t97, -t97 * t109 - t92 * t94, -t99 * t102 + t106, (-t99 * t104 + t112) * t101 + (t99 * t111 + t107) * t100 - pkin(7) * t109 + t95 * pkin(1) + 0; t97 * t108 + t99 * t94, -t94 * t108 + t99 * t97, -t96 * t101, t99 * pkin(7) + qJ(1) + 0 + (pkin(2) * t100 - qJ(3) * t101) * t96; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:00:19
	% EndTime: 2020-11-04 21:00:19
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (54->31), mult. (72->49), div. (0->0), fcn. (93->10), ass. (0->23)
	t116 = cos(pkin(11)) * pkin(3) + pkin(2);
	t120 = sin(pkin(10));
	t134 = t120 * t116;
	t121 = sin(pkin(6));
	t133 = t120 * t121;
	t122 = cos(pkin(10));
	t132 = t121 * t122;
	t125 = sin(qJ(2));
	t131 = t121 * t125;
	t130 = t122 * t116;
	t123 = cos(pkin(6));
	t124 = qJ(3) + pkin(8);
	t129 = t123 * t124;
	t128 = t123 * t125;
	t126 = cos(qJ(2));
	t127 = t123 * t126;
	t119 = pkin(11) + qJ(4);
	t118 = cos(t119);
	t117 = sin(t119);
	t115 = sin(pkin(11)) * pkin(3) + pkin(7);
	t114 = -t120 * t128 + t122 * t126;
	t113 = t120 * t126 + t122 * t128;
	t1 = [t114 * t118 + t117 * t133, -t114 * t117 + t118 * t133, t120 * t127 + t122 * t125, (t120 * t129 + t130) * t126 + (t122 * t124 - t123 * t134) * t125 + t115 * t133 + t122 * pkin(1) + 0; t113 * t118 - t117 * t132, -t113 * t117 - t118 * t132, t120 * t125 - t122 * t127, (-t122 * t129 + t134) * t126 + (t120 * t124 + t123 * t130) * t125 - t115 * t132 + t120 * pkin(1) + 0; t123 * t117 + t118 * t131, -t117 * t131 + t123 * t118, -t121 * t126, t115 * t123 + qJ(1) + 0 + (t116 * t125 - t124 * t126) * t121; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:00:19
	% EndTime: 2020-11-04 21:00:19
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (82->41), mult. (108->55), div. (0->0), fcn. (139->10), ass. (0->29)
	t144 = cos(pkin(11)) * pkin(3) + pkin(2);
	t148 = sin(pkin(10));
	t162 = t148 * t144;
	t149 = sin(pkin(6));
	t161 = t148 * t149;
	t150 = cos(pkin(10));
	t160 = t149 * t150;
	t153 = sin(qJ(2));
	t159 = t149 * t153;
	t158 = t150 * t144;
	t151 = cos(pkin(6));
	t152 = qJ(3) + pkin(8);
	t157 = t151 * t152;
	t156 = t151 * t153;
	t154 = cos(qJ(2));
	t155 = t151 * t154;
	t147 = pkin(11) + qJ(4);
	t146 = cos(t147);
	t145 = sin(t147);
	t143 = sin(pkin(11)) * pkin(3) + pkin(7);
	t142 = -t148 * t156 + t150 * t154;
	t141 = t148 * t154 + t150 * t156;
	t140 = t151 * t145 + t146 * t159;
	t139 = t145 * t159 - t151 * t146;
	t138 = -t142 * t145 + t146 * t161;
	t137 = t142 * t146 + t145 * t161;
	t136 = t141 * t146 - t145 * t160;
	t135 = t141 * t145 + t146 * t160;
	t1 = [t148 * t155 + t150 * t153, -t137, -t138, t137 * pkin(4) - t138 * qJ(5) + (t148 * t157 + t158) * t154 + (t150 * t152 - t151 * t162) * t153 + t143 * t161 + t150 * pkin(1) + 0; t148 * t153 - t150 * t155, -t136, t135, t136 * pkin(4) + t135 * qJ(5) + (-t150 * t157 + t162) * t154 + (t148 * t152 + t151 * t158) * t153 - t143 * t160 + t148 * pkin(1) + 0; -t149 * t154, -t140, t139, t140 * pkin(4) + t139 * qJ(5) + t143 * t151 + qJ(1) + 0 + (t144 * t153 - t152 * t154) * t149; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:00:19
	% EndTime: 2020-11-04 21:00:19
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (119->46), mult. (170->69), div. (0->0), fcn. (221->12), ass. (0->35)
	t196 = pkin(9) + pkin(4);
	t174 = cos(pkin(11)) * pkin(3) + pkin(2);
	t178 = sin(pkin(10));
	t195 = t178 * t174;
	t179 = sin(pkin(6));
	t194 = t178 * t179;
	t180 = cos(pkin(10));
	t193 = t179 * t180;
	t184 = sin(qJ(2));
	t192 = t179 * t184;
	t186 = cos(qJ(2));
	t191 = t179 * t186;
	t190 = t180 * t174;
	t181 = cos(pkin(6));
	t182 = qJ(3) + pkin(8);
	t189 = t181 * t182;
	t188 = t181 * t184;
	t187 = t181 * t186;
	t185 = cos(qJ(6));
	t183 = sin(qJ(6));
	t177 = pkin(11) + qJ(4);
	t176 = cos(t177);
	t175 = sin(t177);
	t173 = sin(pkin(11)) * pkin(3) + pkin(7);
	t172 = -t178 * t188 + t180 * t186;
	t171 = t178 * t187 + t180 * t184;
	t170 = t178 * t186 + t180 * t188;
	t169 = t178 * t184 - t180 * t187;
	t168 = t181 * t175 + t176 * t192;
	t167 = t175 * t192 - t181 * t176;
	t166 = -t172 * t175 + t176 * t194;
	t165 = t172 * t176 + t175 * t194;
	t164 = t170 * t176 - t175 * t193;
	t163 = t170 * t175 + t176 * t193;
	t1 = [-t166 * t183 + t171 * t185, -t166 * t185 - t171 * t183, t165, t171 * pkin(5) - t166 * qJ(5) + (t178 * t189 + t190) * t186 + (t180 * t182 - t181 * t195) * t184 + t173 * t194 + t180 * pkin(1) + 0 + t196 * t165; t163 * t183 + t169 * t185, t163 * t185 - t169 * t183, t164, t169 * pkin(5) + t163 * qJ(5) + (-t180 * t189 + t195) * t186 + (t178 * t182 + t181 * t190) * t184 - t173 * t193 + t178 * pkin(1) + 0 + t196 * t164; t167 * t183 - t185 * t191, t167 * t185 + t183 * t191, t168, t167 * qJ(5) + t173 * t181 + qJ(1) + 0 + t196 * t168 + (t174 * t184 + (-pkin(5) - t182) * t186) * t179; 0, 0, 0, 1;];
	Tc_mdh = t1;
end