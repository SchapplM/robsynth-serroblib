% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRPRRR3 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:04
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PRPRRR3_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR3_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRR3_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR3_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:04:28
	% EndTime: 2020-11-04 21:04:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:04:28
	% EndTime: 2020-11-04 21:04:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t86 = cos(pkin(11));
	t85 = sin(pkin(11));
	t1 = [t86, -t85, 0, 0; t85, t86, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:04:28
	% EndTime: 2020-11-04 21:04:28
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t87 = sin(pkin(11));
	t88 = sin(pkin(6));
	t96 = t87 * t88;
	t89 = cos(pkin(11));
	t95 = t89 * t88;
	t90 = cos(pkin(6));
	t91 = sin(qJ(2));
	t94 = t90 * t91;
	t92 = cos(qJ(2));
	t93 = t90 * t92;
	t1 = [-t87 * t94 + t89 * t92, -t87 * t93 - t89 * t91, t96, t89 * pkin(1) + pkin(7) * t96 + 0; t87 * t92 + t89 * t94, -t87 * t91 + t89 * t93, -t95, t87 * pkin(1) - pkin(7) * t95 + 0; t88 * t91, t88 * t92, t90, t90 * pkin(7) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:04:28
	% EndTime: 2020-11-04 21:04:28
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (29->27), mult. (64->48), div. (0->0), fcn. (85->8), ass. (0->18)
	t100 = sin(pkin(11));
	t101 = sin(pkin(6));
	t113 = t100 * t101;
	t104 = cos(pkin(6));
	t112 = t100 * t104;
	t103 = cos(pkin(11));
	t111 = t101 * t103;
	t105 = sin(qJ(2));
	t110 = t101 * t105;
	t109 = t103 * t104;
	t108 = t104 * t105;
	t106 = cos(qJ(2));
	t107 = t104 * t106;
	t102 = cos(pkin(12));
	t99 = sin(pkin(12));
	t98 = t100 * t106 + t103 * t108;
	t97 = t100 * t108 - t103 * t106;
	t1 = [-t97 * t102 + t99 * t113, t102 * t113 + t97 * t99, t100 * t107 + t103 * t105, (t103 * pkin(2) + qJ(3) * t112) * t106 + (-pkin(2) * t112 + t103 * qJ(3)) * t105 + pkin(7) * t113 + t103 * pkin(1) + 0; t98 * t102 - t99 * t111, -t102 * t111 - t98 * t99, t100 * t105 - t103 * t107, (t100 * pkin(2) - qJ(3) * t109) * t106 + (pkin(2) * t109 + t100 * qJ(3)) * t105 - pkin(7) * t111 + t100 * pkin(1) + 0; t102 * t110 + t104 * t99, t104 * t102 - t99 * t110, -t101 * t106, t104 * pkin(7) + qJ(1) + 0 + (pkin(2) * t105 - qJ(3) * t106) * t101; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:04:28
	% EndTime: 2020-11-04 21:04:28
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (54->31), mult. (72->49), div. (0->0), fcn. (93->10), ass. (0->23)
	t117 = cos(pkin(12)) * pkin(3) + pkin(2);
	t121 = sin(pkin(11));
	t135 = t121 * t117;
	t122 = sin(pkin(6));
	t134 = t121 * t122;
	t123 = cos(pkin(11));
	t133 = t122 * t123;
	t126 = sin(qJ(2));
	t132 = t122 * t126;
	t131 = t123 * t117;
	t124 = cos(pkin(6));
	t125 = qJ(3) + pkin(8);
	t130 = t124 * t125;
	t129 = t124 * t126;
	t127 = cos(qJ(2));
	t128 = t124 * t127;
	t120 = pkin(12) + qJ(4);
	t119 = cos(t120);
	t118 = sin(t120);
	t116 = sin(pkin(12)) * pkin(3) + pkin(7);
	t115 = t121 * t127 + t123 * t129;
	t114 = t121 * t129 - t123 * t127;
	t1 = [-t114 * t119 + t118 * t134, t114 * t118 + t119 * t134, t121 * t128 + t123 * t126, (t121 * t130 + t131) * t127 + (t123 * t125 - t124 * t135) * t126 + t116 * t134 + t123 * pkin(1) + 0; t115 * t119 - t118 * t133, -t115 * t118 - t119 * t133, t121 * t126 - t123 * t128, (-t123 * t130 + t135) * t127 + (t121 * t125 + t124 * t131) * t126 - t116 * t133 + t121 * pkin(1) + 0; t124 * t118 + t119 * t132, -t118 * t132 + t124 * t119, -t122 * t127, t116 * t124 + qJ(1) + 0 + (t117 * t126 - t125 * t127) * t122; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:04:28
	% EndTime: 2020-11-04 21:04:28
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (77->30), mult. (81->42), div. (0->0), fcn. (106->12), ass. (0->23)
	t146 = pkin(12) + qJ(4);
	t158 = pkin(7) + pkin(4) * sin(t146) + sin(pkin(12)) * pkin(3);
	t147 = sin(pkin(11));
	t148 = sin(pkin(6));
	t157 = t147 * t148;
	t149 = cos(pkin(11));
	t156 = t148 * t149;
	t151 = sin(qJ(2));
	t155 = t148 * t151;
	t150 = cos(pkin(6));
	t154 = t150 * t151;
	t152 = cos(qJ(2));
	t153 = t150 * t152;
	t145 = -pkin(9) - pkin(8) - qJ(3);
	t144 = qJ(5) + t146;
	t143 = cos(t144);
	t142 = sin(t144);
	t140 = pkin(4) * cos(t146) + cos(pkin(12)) * pkin(3) + pkin(2);
	t139 = -t147 * t154 + t149 * t152;
	t138 = t147 * t153 + t149 * t151;
	t137 = t147 * t152 + t149 * t154;
	t136 = t147 * t151 - t149 * t153;
	t1 = [t139 * t143 + t142 * t157, -t139 * t142 + t143 * t157, t138, t149 * pkin(1) - t138 * t145 + t139 * t140 + t158 * t157 + 0; t137 * t143 - t142 * t156, -t137 * t142 - t143 * t156, t136, t147 * pkin(1) - t136 * t145 + t137 * t140 - t158 * t156 + 0; t150 * t142 + t143 * t155, -t142 * t155 + t150 * t143, -t148 * t152, qJ(1) + 0 + t158 * t150 + (t140 * t151 + t145 * t152) * t148; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:04:28
	% EndTime: 2020-11-04 21:04:28
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (141->41), mult. (151->60), div. (0->0), fcn. (199->14), ass. (0->32)
	t175 = pkin(12) + qJ(4);
	t190 = pkin(7) + pkin(4) * sin(t175) + sin(pkin(12)) * pkin(3);
	t176 = sin(pkin(11));
	t177 = sin(pkin(6));
	t189 = t176 * t177;
	t178 = cos(pkin(11));
	t188 = t177 * t178;
	t181 = sin(qJ(2));
	t187 = t177 * t181;
	t183 = cos(qJ(2));
	t186 = t177 * t183;
	t179 = cos(pkin(6));
	t185 = t179 * t181;
	t184 = t179 * t183;
	t182 = cos(qJ(6));
	t180 = sin(qJ(6));
	t174 = -pkin(9) - pkin(8) - qJ(3);
	t173 = qJ(5) + t175;
	t172 = cos(t173);
	t171 = sin(t173);
	t169 = pkin(4) * cos(t175) + cos(pkin(12)) * pkin(3) + pkin(2);
	t168 = -t176 * t185 + t178 * t183;
	t167 = t176 * t184 + t178 * t181;
	t166 = t176 * t183 + t178 * t185;
	t165 = t176 * t181 - t178 * t184;
	t164 = t179 * t171 + t172 * t187;
	t163 = t171 * t187 - t179 * t172;
	t162 = t168 * t172 + t171 * t189;
	t161 = t168 * t171 - t172 * t189;
	t160 = t166 * t172 - t171 * t188;
	t159 = t166 * t171 + t172 * t188;
	t1 = [t162 * t182 + t167 * t180, -t162 * t180 + t167 * t182, t161, t178 * pkin(1) + t162 * pkin(5) + t161 * pkin(10) - t167 * t174 + t168 * t169 + t190 * t189 + 0; t160 * t182 + t165 * t180, -t160 * t180 + t165 * t182, t159, t176 * pkin(1) + t160 * pkin(5) + t159 * pkin(10) - t165 * t174 + t166 * t169 - t190 * t188 + 0; t164 * t182 - t180 * t186, -t164 * t180 - t182 * t186, t163, t164 * pkin(5) + t163 * pkin(10) + qJ(1) + 0 + t190 * t179 + (t169 * t181 + t174 * t183) * t177; 0, 0, 0, 1;];
	Tc_mdh = t1;
end