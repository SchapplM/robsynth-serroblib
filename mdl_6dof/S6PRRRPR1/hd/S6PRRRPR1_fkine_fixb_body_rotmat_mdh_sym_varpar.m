% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRRPR1 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:14
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PRRRPR1_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR1_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPR1_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR1_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:14:38
	% EndTime: 2020-11-04 21:14:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:14:38
	% EndTime: 2020-11-04 21:14:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t84 = cos(pkin(11));
	t83 = sin(pkin(11));
	t1 = [t84, -t83, 0, 0; t83, t84, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:14:38
	% EndTime: 2020-11-04 21:14:38
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t85 = sin(pkin(11));
	t86 = sin(pkin(6));
	t94 = t85 * t86;
	t87 = cos(pkin(11));
	t93 = t87 * t86;
	t88 = cos(pkin(6));
	t89 = sin(qJ(2));
	t92 = t88 * t89;
	t90 = cos(qJ(2));
	t91 = t88 * t90;
	t1 = [-t85 * t92 + t87 * t90, -t85 * t91 - t87 * t89, t94, t87 * pkin(1) + pkin(7) * t94 + 0; t85 * t90 + t87 * t92, -t85 * t89 + t87 * t91, -t93, t85 * pkin(1) - pkin(7) * t93 + 0; t86 * t89, t86 * t90, t88, t88 * pkin(7) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:14:38
	% EndTime: 2020-11-04 21:14:38
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (29->27), mult. (64->48), div. (0->0), fcn. (85->8), ass. (0->18)
	t98 = sin(pkin(6));
	t111 = t98 * pkin(7);
	t100 = cos(pkin(6));
	t97 = sin(pkin(11));
	t110 = t100 * t97;
	t99 = cos(pkin(11));
	t109 = t100 * t99;
	t101 = sin(qJ(3));
	t108 = t101 * t98;
	t103 = cos(qJ(3));
	t107 = t103 * t98;
	t102 = sin(qJ(2));
	t106 = t100 * t102;
	t104 = cos(qJ(2));
	t105 = t100 * t104;
	t96 = t97 * t104 + t99 * t106;
	t95 = -t99 * t104 + t97 * t106;
	t1 = [-t95 * t103 + t97 * t108, t95 * t101 + t97 * t107, t99 * t102 + t97 * t105, (t99 * pkin(2) + pkin(8) * t110) * t104 + (-pkin(2) * t110 + t99 * pkin(8)) * t102 + t97 * t111 + t99 * pkin(1) + 0; t96 * t103 - t99 * t108, -t96 * t101 - t99 * t107, t97 * t102 - t99 * t105, (t97 * pkin(2) - pkin(8) * t109) * t104 + (pkin(2) * t109 + t97 * pkin(8)) * t102 - t99 * t111 + t97 * pkin(1) + 0; t100 * t101 + t102 * t107, t100 * t103 - t102 * t108, -t98 * t104, t100 * pkin(7) + qJ(1) + 0 + (pkin(2) * t102 - pkin(8) * t104) * t98; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:14:38
	% EndTime: 2020-11-04 21:14:38
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (50->27), mult. (75->40), div. (0->0), fcn. (100->10), ass. (0->22)
	t121 = sin(pkin(11));
	t122 = sin(pkin(6));
	t134 = t121 * t122;
	t123 = cos(pkin(11));
	t133 = t122 * t123;
	t126 = sin(qJ(2));
	t132 = t122 * t126;
	t124 = cos(pkin(6));
	t131 = t124 * t126;
	t127 = cos(qJ(2));
	t130 = t124 * t127;
	t129 = pkin(3) * sin(qJ(3)) + pkin(7);
	t128 = pkin(9) + pkin(8);
	t120 = qJ(3) + qJ(4);
	t119 = cos(t120);
	t118 = sin(t120);
	t117 = cos(qJ(3)) * pkin(3) + pkin(2);
	t115 = t121 * t130 + t123 * t126;
	t114 = t121 * t127 + t123 * t131;
	t113 = t121 * t126 - t123 * t130;
	t112 = t121 * t131 - t123 * t127;
	t1 = [-t112 * t119 + t118 * t134, t112 * t118 + t119 * t134, t115, t123 * pkin(1) - t112 * t117 + t115 * t128 + t129 * t134 + 0; t114 * t119 - t118 * t133, -t114 * t118 - t119 * t133, t113, t121 * pkin(1) + t113 * t128 + t114 * t117 - t129 * t133 + 0; t124 * t118 + t119 * t132, -t118 * t132 + t124 * t119, -t122 * t127, qJ(1) + 0 + t129 * t124 + (t117 * t126 - t127 * t128) * t122; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:14:38
	% EndTime: 2020-11-04 21:14:38
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (77->30), mult. (81->42), div. (0->0), fcn. (106->12), ass. (0->23)
	t145 = qJ(3) + qJ(4);
	t157 = pkin(7) + pkin(4) * sin(t145) + pkin(3) * sin(qJ(3));
	t146 = sin(pkin(11));
	t147 = sin(pkin(6));
	t156 = t146 * t147;
	t148 = cos(pkin(11));
	t155 = t147 * t148;
	t150 = sin(qJ(2));
	t154 = t147 * t150;
	t149 = cos(pkin(6));
	t153 = t149 * t150;
	t151 = cos(qJ(2));
	t152 = t149 * t151;
	t144 = -qJ(5) - pkin(9) - pkin(8);
	t143 = pkin(12) + t145;
	t142 = cos(t143);
	t141 = sin(t143);
	t139 = pkin(4) * cos(t145) + cos(qJ(3)) * pkin(3) + pkin(2);
	t138 = -t146 * t153 + t148 * t151;
	t137 = t146 * t152 + t148 * t150;
	t136 = t146 * t151 + t148 * t153;
	t135 = t146 * t150 - t148 * t152;
	t1 = [t138 * t142 + t141 * t156, -t138 * t141 + t142 * t156, t137, t148 * pkin(1) - t137 * t144 + t138 * t139 + t157 * t156 + 0; t136 * t142 - t141 * t155, -t136 * t141 - t142 * t155, t135, t146 * pkin(1) - t135 * t144 + t136 * t139 - t157 * t155 + 0; t149 * t141 + t142 * t154, -t141 * t154 + t149 * t142, -t147 * t151, qJ(1) + 0 + t157 * t149 + (t139 * t150 + t144 * t151) * t147; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:14:38
	% EndTime: 2020-11-04 21:14:38
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (141->41), mult. (151->60), div. (0->0), fcn. (199->14), ass. (0->32)
	t174 = qJ(3) + qJ(4);
	t189 = pkin(7) + pkin(4) * sin(t174) + pkin(3) * sin(qJ(3));
	t175 = sin(pkin(11));
	t176 = sin(pkin(6));
	t188 = t175 * t176;
	t177 = cos(pkin(11));
	t187 = t176 * t177;
	t180 = sin(qJ(2));
	t186 = t176 * t180;
	t182 = cos(qJ(2));
	t185 = t176 * t182;
	t178 = cos(pkin(6));
	t184 = t178 * t180;
	t183 = t178 * t182;
	t181 = cos(qJ(6));
	t179 = sin(qJ(6));
	t173 = -qJ(5) - pkin(9) - pkin(8);
	t172 = pkin(12) + t174;
	t171 = cos(t172);
	t170 = sin(t172);
	t168 = pkin(4) * cos(t174) + cos(qJ(3)) * pkin(3) + pkin(2);
	t167 = -t175 * t184 + t177 * t182;
	t166 = t175 * t183 + t177 * t180;
	t165 = t175 * t182 + t177 * t184;
	t164 = t175 * t180 - t177 * t183;
	t163 = t178 * t170 + t171 * t186;
	t162 = t170 * t186 - t178 * t171;
	t161 = t167 * t171 + t170 * t188;
	t160 = t167 * t170 - t171 * t188;
	t159 = t165 * t171 - t170 * t187;
	t158 = t165 * t170 + t171 * t187;
	t1 = [t161 * t181 + t166 * t179, -t161 * t179 + t166 * t181, t160, t177 * pkin(1) + t161 * pkin(5) + t160 * pkin(10) - t166 * t173 + t167 * t168 + t189 * t188 + 0; t159 * t181 + t164 * t179, -t159 * t179 + t164 * t181, t158, t175 * pkin(1) + t159 * pkin(5) + t158 * pkin(10) - t164 * t173 + t165 * t168 - t189 * t187 + 0; t163 * t181 - t179 * t185, -t163 * t179 - t181 * t185, t162, t163 * pkin(5) + t162 * pkin(10) + qJ(1) + 0 + t189 * t178 + (t168 * t180 + t173 * t182) * t176; 0, 0, 0, 1;];
	Tc_mdh = t1;
end