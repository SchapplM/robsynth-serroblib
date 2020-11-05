% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRPRR4 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:00
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5PRPRR4_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR4_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPRR4_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR4_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:00:15
	% EndTime: 2020-11-04 20:00:15
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:00:15
	% EndTime: 2020-11-04 20:00:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t89 = cos(pkin(9));
	t88 = sin(pkin(9));
	t1 = [t89, -t88, 0, 0; t88, t89, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:00:15
	% EndTime: 2020-11-04 20:00:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t90 = sin(pkin(9));
	t91 = sin(pkin(5));
	t99 = t90 * t91;
	t92 = cos(pkin(9));
	t98 = t92 * t91;
	t93 = cos(pkin(5));
	t94 = sin(qJ(2));
	t97 = t93 * t94;
	t95 = cos(qJ(2));
	t96 = t93 * t95;
	t1 = [-t90 * t97 + t92 * t95, -t90 * t96 - t92 * t94, t99, t92 * pkin(1) + pkin(6) * t99 + 0; t90 * t95 + t92 * t97, -t90 * t94 + t92 * t96, -t98, t90 * pkin(1) - pkin(6) * t98 + 0; t91 * t94, t91 * t95, t93, t93 * pkin(6) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:00:15
	% EndTime: 2020-11-04 20:00:15
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (53->22), mult. (37->26), div. (0->0), fcn. (44->12), ass. (0->20)
	t119 = pkin(2) * sin(qJ(2));
	t112 = qJ(2) + pkin(10);
	t117 = pkin(6) + qJ(3);
	t116 = cos(pkin(5));
	t115 = cos(pkin(9));
	t114 = sin(pkin(5));
	t113 = sin(pkin(9));
	t111 = cos(t112);
	t110 = sin(t112);
	t109 = pkin(5) - t112;
	t108 = pkin(5) + t112;
	t107 = cos(qJ(2)) * pkin(2) + pkin(1);
	t106 = cos(t108);
	t105 = sin(t109);
	t104 = cos(t109) / 0.2e1;
	t103 = sin(t108) / 0.2e1;
	t102 = -t114 * t117 + t116 * t119;
	t101 = t106 / 0.2e1 + t104;
	t100 = t103 - t105 / 0.2e1;
	t1 = [-t113 * t100 + t115 * t111, -t113 * t101 - t115 * t110, t113 * t114, -t113 * t102 + t115 * t107 + 0; t115 * t100 + t113 * t111, t115 * t101 - t113 * t110, -t115 * t114, t115 * t102 + t113 * t107 + 0; t104 - t106 / 0.2e1, t105 / 0.2e1 + t103, t116, t114 * t119 + t116 * t117 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:00:15
	% EndTime: 2020-11-04 20:00:15
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (73->35), mult. (85->53), div. (0->0), fcn. (105->16), ass. (0->27)
	t131 = sin(pkin(9));
	t135 = cos(pkin(5));
	t145 = t131 * t135;
	t132 = sin(pkin(5));
	t136 = pkin(6) + qJ(3);
	t144 = t132 * t136;
	t137 = sin(qJ(4));
	t143 = t132 * t137;
	t139 = cos(qJ(4));
	t142 = t132 * t139;
	t134 = cos(pkin(9));
	t141 = t134 * t135;
	t129 = qJ(2) + pkin(10);
	t140 = cos(qJ(2));
	t138 = sin(qJ(2));
	t133 = cos(pkin(10));
	t130 = sin(pkin(10));
	t128 = cos(t129);
	t127 = sin(t129);
	t126 = pkin(5) - t129;
	t125 = pkin(5) + t129;
	t124 = -t130 * pkin(3) + t133 * pkin(7);
	t123 = t133 * pkin(3) + t130 * pkin(7) + pkin(2);
	t122 = cos(t125) + cos(t126);
	t121 = t127 * t141 + t131 * t128;
	t120 = t127 * t145 - t134 * t128;
	t1 = [-t120 * t139 + t131 * t143, t120 * t137 + t131 * t142, t134 * t127 + t131 * t122 / 0.2e1, (t134 * t123 + t124 * t145) * t140 + (-t123 * t145 + t134 * t124) * t138 + t131 * t144 + t134 * pkin(1) + 0; t121 * t139 - t134 * t143, -t121 * t137 - t134 * t142, t131 * t127 - t134 * t122 / 0.2e1, (t131 * t123 - t124 * t141) * t140 + (t123 * t141 + t131 * t124) * t138 - t134 * t144 + t131 * pkin(1) + 0; t127 * t142 + t135 * t137, -t127 * t143 + t135 * t139, -sin(t126) / 0.2e1 - sin(t125) / 0.2e1, t135 * t136 + qJ(1) + 0 + (t123 * t138 - t124 * t140) * t132; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:00:15
	% EndTime: 2020-11-04 20:00:15
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (134->48), mult. (156->71), div. (0->0), fcn. (198->18), ass. (0->38)
	t166 = sin(pkin(9));
	t170 = cos(pkin(5));
	t182 = t166 * t170;
	t167 = sin(pkin(5));
	t171 = pkin(6) + qJ(3);
	t181 = t167 * t171;
	t173 = sin(qJ(4));
	t180 = t167 * t173;
	t176 = cos(qJ(4));
	t179 = t167 * t176;
	t169 = cos(pkin(9));
	t178 = t169 * t170;
	t164 = qJ(2) + pkin(10);
	t177 = cos(qJ(2));
	t175 = cos(qJ(5));
	t174 = sin(qJ(2));
	t172 = sin(qJ(5));
	t168 = cos(pkin(10));
	t165 = sin(pkin(10));
	t163 = cos(t164);
	t162 = sin(t164);
	t161 = pkin(5) - t164;
	t160 = pkin(5) + t164;
	t159 = -t165 * pkin(3) + t168 * pkin(7);
	t158 = t168 * pkin(3) + t165 * pkin(7) + pkin(2);
	t157 = cos(t160) + cos(t161);
	t156 = -sin(t161) / 0.2e1 - sin(t160) / 0.2e1;
	t155 = t162 * t179 + t170 * t173;
	t154 = t162 * t180 - t170 * t176;
	t153 = t162 * t178 + t166 * t163;
	t152 = t162 * t182 - t169 * t163;
	t151 = t169 * t162 + t166 * t157 / 0.2e1;
	t150 = t166 * t162 - t169 * t157 / 0.2e1;
	t149 = -t152 * t176 + t166 * t180;
	t148 = t153 * t176 - t169 * t180;
	t147 = t153 * t173 + t169 * t179;
	t146 = t152 * t173 + t166 * t179;
	t1 = [t149 * t175 + t151 * t172, -t149 * t172 + t151 * t175, -t146, t149 * pkin(4) - t146 * pkin(8) + (t169 * t158 + t159 * t182) * t177 + (-t158 * t182 + t169 * t159) * t174 + t166 * t181 + t169 * pkin(1) + 0; t148 * t175 + t150 * t172, -t148 * t172 + t150 * t175, t147, t148 * pkin(4) + t147 * pkin(8) + (t166 * t158 - t159 * t178) * t177 + (t158 * t178 + t166 * t159) * t174 - t169 * t181 + t166 * pkin(1) + 0; t155 * t175 + t156 * t172, -t155 * t172 + t156 * t175, t154, t155 * pkin(4) + t154 * pkin(8) + t170 * t171 + qJ(1) + 0 + (t158 * t174 - t159 * t177) * t167; 0, 0, 0, 1;];
	Tc_mdh = t1;
end