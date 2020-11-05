% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRPRRR5 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:05
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PRPRRR5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRR5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:05:13
	% EndTime: 2020-11-04 21:05:13
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:05:13
	% EndTime: 2020-11-04 21:05:13
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t85 = cos(pkin(11));
	t84 = sin(pkin(11));
	t1 = [t85, -t84, 0, 0; t84, t85, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:05:13
	% EndTime: 2020-11-04 21:05:13
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t86 = sin(pkin(11));
	t87 = sin(pkin(6));
	t95 = t86 * t87;
	t88 = cos(pkin(11));
	t94 = t88 * t87;
	t89 = cos(pkin(6));
	t90 = sin(qJ(2));
	t93 = t89 * t90;
	t91 = cos(qJ(2));
	t92 = t89 * t91;
	t1 = [-t86 * t93 + t88 * t91, -t86 * t92 - t88 * t90, t95, t88 * pkin(1) + pkin(7) * t95 + 0; t86 * t91 + t88 * t93, -t86 * t90 + t88 * t92, -t94, t86 * pkin(1) - pkin(7) * t94 + 0; t87 * t90, t87 * t91, t89, t89 * pkin(7) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:05:13
	% EndTime: 2020-11-04 21:05:13
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (23->23), mult. (43->32), div. (0->0), fcn. (56->6), ass. (0->17)
	t96 = sin(pkin(11));
	t111 = t96 * pkin(2);
	t98 = cos(pkin(11));
	t110 = t98 * pkin(2);
	t97 = sin(pkin(6));
	t109 = t96 * t97;
	t108 = t98 * t97;
	t107 = t96 * qJ(3);
	t100 = sin(qJ(2));
	t106 = t96 * t100;
	t101 = cos(qJ(2));
	t105 = t96 * t101;
	t104 = t98 * qJ(3);
	t103 = t98 * t100;
	t102 = t98 * t101;
	t99 = cos(pkin(6));
	t1 = [t109, t99 * t106 - t102, t99 * t105 + t103, (t99 * t107 + t110) * t101 + (-t99 * t111 + t104) * t100 + pkin(7) * t109 + t98 * pkin(1) + 0; -t108, -t99 * t103 - t105, -t99 * t102 + t106, (-t99 * t104 + t111) * t101 + (t99 * t110 + t107) * t100 - pkin(7) * t108 + t96 * pkin(1) + 0; t99, -t97 * t100, -t97 * t101, t99 * pkin(7) + qJ(1) + 0 + (pkin(2) * t100 - qJ(3) * t101) * t97; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:05:13
	% EndTime: 2020-11-04 21:05:13
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (36->28), mult. (64->48), div. (0->0), fcn. (85->8), ass. (0->22)
	t114 = sin(pkin(11));
	t132 = t114 * qJ(3);
	t115 = sin(pkin(6));
	t118 = sin(qJ(4));
	t131 = t115 * t118;
	t120 = cos(qJ(4));
	t130 = t115 * t120;
	t121 = cos(qJ(2));
	t129 = t115 * t121;
	t122 = pkin(3) + pkin(7);
	t128 = t115 * t122;
	t116 = cos(pkin(11));
	t127 = t116 * qJ(3);
	t117 = cos(pkin(6));
	t119 = sin(qJ(2));
	t126 = t117 * t119;
	t125 = t117 * t121;
	t123 = pkin(2) + pkin(8);
	t124 = t117 * t123;
	t113 = t114 * t125 + t116 * t119;
	t112 = t114 * t119 - t116 * t125;
	t1 = [t113 * t118 + t114 * t130, t113 * t120 - t114 * t131, -t114 * t126 + t116 * t121, (t116 * t123 + t117 * t132) * t121 + (-t114 * t124 + t127) * t119 + t114 * t128 + t116 * pkin(1) + 0; t112 * t118 - t116 * t130, t112 * t120 + t116 * t131, t114 * t121 + t116 * t126, (t114 * t123 - t117 * t127) * t121 + (t116 * t124 + t132) * t119 - t116 * t128 + t114 * pkin(1) + 0; t117 * t120 - t118 * t129, -t117 * t118 - t120 * t129, t115 * t119, t122 * t117 + qJ(1) + 0 + (-qJ(3) * t121 + t119 * t123) * t115; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:05:13
	% EndTime: 2020-11-04 21:05:13
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (61->26), mult. (83->40), div. (0->0), fcn. (104->10), ass. (0->23)
	t140 = sin(pkin(6));
	t142 = cos(pkin(6));
	t137 = pkin(2) + pkin(8) + pkin(9);
	t144 = sin(qJ(2));
	t146 = cos(qJ(2));
	t151 = sin(qJ(4)) * pkin(4) + qJ(3);
	t149 = -t137 * t144 + t151 * t146;
	t150 = cos(qJ(4)) * pkin(4) + pkin(3) + pkin(7);
	t160 = t150 * t140 + t149 * t142;
	t139 = sin(pkin(11));
	t156 = t139 * t140;
	t141 = cos(pkin(11));
	t155 = t140 * t141;
	t154 = t140 * t146;
	t153 = t142 * t144;
	t152 = t142 * t146;
	t148 = t137 * t146 + t151 * t144 + pkin(1);
	t138 = qJ(4) + qJ(5);
	t136 = cos(t138);
	t135 = sin(t138);
	t134 = t139 * t152 + t141 * t144;
	t133 = t139 * t144 - t141 * t152;
	t1 = [t134 * t135 + t136 * t156, t134 * t136 - t135 * t156, -t139 * t153 + t141 * t146, t160 * t139 + t148 * t141 + 0; t133 * t135 - t136 * t155, t133 * t136 + t135 * t155, t139 * t146 + t141 * t153, t148 * t139 - t160 * t141 + 0; -t135 * t154 + t142 * t136, -t142 * t135 - t136 * t154, t140 * t144, -t149 * t140 + t150 * t142 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:05:13
	% EndTime: 2020-11-04 21:05:13
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (108->39), mult. (153->58), div. (0->0), fcn. (197->12), ass. (0->34)
	t176 = sin(pkin(6));
	t178 = cos(pkin(6));
	t173 = pkin(2) + pkin(8) + pkin(9);
	t181 = sin(qJ(2));
	t184 = cos(qJ(2));
	t189 = sin(qJ(4)) * pkin(4) + qJ(3);
	t187 = -t173 * t181 + t189 * t184;
	t188 = cos(qJ(4)) * pkin(4) + pkin(3) + pkin(7);
	t199 = t188 * t176 + t187 * t178;
	t175 = sin(pkin(11));
	t195 = t175 * t176;
	t177 = cos(pkin(11));
	t194 = t176 * t177;
	t193 = t176 * t181;
	t192 = t176 * t184;
	t191 = t178 * t181;
	t190 = t178 * t184;
	t186 = t173 * t184 + t189 * t181 + pkin(1);
	t182 = cos(qJ(6));
	t179 = sin(qJ(6));
	t174 = qJ(4) + qJ(5);
	t172 = cos(t174);
	t171 = sin(t174);
	t170 = -t175 * t191 + t177 * t184;
	t169 = t175 * t190 + t177 * t181;
	t168 = t175 * t184 + t177 * t191;
	t167 = t175 * t181 - t177 * t190;
	t166 = -t171 * t192 + t178 * t172;
	t165 = t178 * t171 + t172 * t192;
	t164 = t167 * t171 - t172 * t194;
	t163 = t167 * t172 + t171 * t194;
	t162 = t169 * t171 + t172 * t195;
	t161 = -t169 * t172 + t171 * t195;
	t1 = [t162 * t182 + t170 * t179, -t162 * t179 + t170 * t182, t161, t162 * pkin(5) + t161 * pkin(10) + t199 * t175 + t186 * t177 + 0; t164 * t182 + t168 * t179, -t164 * t179 + t168 * t182, -t163, t164 * pkin(5) - t163 * pkin(10) + t186 * t175 - t199 * t177 + 0; t166 * t182 + t179 * t193, -t166 * t179 + t182 * t193, t165, t166 * pkin(5) + t165 * pkin(10) - t187 * t176 + t188 * t178 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end