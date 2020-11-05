% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPRR10 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:37
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RRPRR10_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR10_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRR10_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR10_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:37:56
	% EndTime: 2020-11-04 20:37:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:37:56
	% EndTime: 2020-11-04 20:37:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t95 = cos(qJ(1));
	t94 = sin(qJ(1));
	t1 = [t95, -t94, 0, 0; t94, t95, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:37:56
	% EndTime: 2020-11-04 20:37:56
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t96 = sin(pkin(5));
	t99 = sin(qJ(1));
	t107 = t99 * t96;
	t98 = sin(qJ(2));
	t106 = t99 * t98;
	t101 = cos(qJ(1));
	t105 = t101 * t96;
	t104 = t101 * t98;
	t100 = cos(qJ(2));
	t103 = t99 * t100;
	t102 = t101 * t100;
	t97 = cos(pkin(5));
	t1 = [-t97 * t106 + t102, -t97 * t103 - t104, t107, t101 * pkin(1) + pkin(7) * t107 + 0; t97 * t104 + t103, t97 * t102 - t106, -t105, t99 * pkin(1) - pkin(7) * t105 + 0; t96 * t98, t96 * t100, t97, t97 * pkin(7) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:37:56
	% EndTime: 2020-11-04 20:37:56
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (53->22), mult. (33->28), div. (0->0), fcn. (44->12), ass. (0->20)
	t127 = pkin(2) * sin(qJ(2));
	t120 = qJ(2) + pkin(10);
	t126 = cos(qJ(1));
	t125 = sin(qJ(1));
	t123 = pkin(7) + qJ(3);
	t122 = cos(pkin(5));
	t121 = sin(pkin(5));
	t119 = cos(t120);
	t118 = sin(t120);
	t117 = pkin(5) - t120;
	t116 = pkin(5) + t120;
	t115 = cos(qJ(2)) * pkin(2) + pkin(1);
	t114 = cos(t117);
	t113 = cos(t116);
	t112 = sin(t117);
	t111 = sin(t116);
	t110 = t113 + t114;
	t109 = -t111 + t112;
	t108 = -t121 * t123 + t122 * t127;
	t1 = [t126 * t119 + t125 * t109 / 0.2e1, -t126 * t118 - t125 * t110 / 0.2e1, t125 * t121, -t108 * t125 + t126 * t115 + 0; t125 * t119 - t126 * t109 / 0.2e1, -t125 * t118 + t126 * t110 / 0.2e1, -t126 * t121, t126 * t108 + t125 * t115 + 0; t114 / 0.2e1 - t113 / 0.2e1, t112 / 0.2e1 + t111 / 0.2e1, t122, t121 * t127 + t122 * t123 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:37:56
	% EndTime: 2020-11-04 20:37:56
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (73->30), mult. (81->43), div. (0->0), fcn. (101->16), ass. (0->30)
	t141 = sin(pkin(5));
	t145 = sin(qJ(4));
	t156 = t141 * t145;
	t148 = cos(qJ(4));
	t155 = t141 * t148;
	t150 = cos(qJ(1));
	t154 = t141 * t150;
	t139 = qJ(2) + pkin(10);
	t137 = sin(t139);
	t147 = sin(qJ(1));
	t153 = t147 * t137;
	t152 = t150 * t137;
	t140 = sin(pkin(10));
	t142 = cos(pkin(10));
	t133 = t142 * pkin(3) + t140 * pkin(8) + pkin(2);
	t134 = -t140 * pkin(3) + t142 * pkin(8);
	t146 = sin(qJ(2));
	t149 = cos(qJ(2));
	t151 = t133 * t146 - t134 * t149;
	t144 = pkin(7) + qJ(3);
	t143 = cos(pkin(5));
	t138 = cos(t139);
	t136 = pkin(5) - t139;
	t135 = pkin(5) + t139;
	t132 = cos(t135) + cos(t136);
	t131 = t147 * t138 + t143 * t152;
	t130 = -t150 * t138 + t143 * t153;
	t129 = t133 * t149 + t134 * t146 + pkin(1);
	t128 = t141 * t144 - t151 * t143;
	t1 = [-t130 * t148 + t147 * t156, t130 * t145 + t147 * t155, t152 + t147 * t132 / 0.2e1, t128 * t147 + t129 * t150 + 0; t131 * t148 - t145 * t154, -t131 * t145 - t148 * t154, t153 - t150 * t132 / 0.2e1, -t128 * t150 + t129 * t147 + 0; t137 * t155 + t143 * t145, -t137 * t156 + t143 * t148, -sin(t136) / 0.2e1 - sin(t135) / 0.2e1, t151 * t141 + t143 * t144 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:37:56
	% EndTime: 2020-11-04 20:37:56
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (102->44), mult. (170->73), div. (0->0), fcn. (204->14), ass. (0->39)
	t164 = qJ(2) + pkin(10);
	t162 = sin(t164);
	t168 = cos(pkin(5));
	t194 = t162 * t168;
	t166 = sin(pkin(5));
	t170 = sin(qJ(4));
	t193 = t166 * t170;
	t173 = cos(qJ(5));
	t192 = t166 * t173;
	t169 = sin(qJ(5));
	t172 = sin(qJ(1));
	t191 = t169 * t172;
	t176 = cos(qJ(1));
	t190 = t169 * t176;
	t174 = cos(qJ(4));
	t189 = t172 * t174;
	t188 = t173 * t172;
	t187 = t173 * t176;
	t186 = t174 * t176;
	t185 = t169 * t193;
	t184 = t170 * t192;
	t183 = t174 * t188;
	t182 = t169 * t189;
	t181 = t169 * t186;
	t180 = t173 * t186;
	t165 = sin(pkin(10));
	t167 = cos(pkin(10));
	t178 = pkin(4) * t174 + pkin(9) * t170 + pkin(3);
	t159 = t165 * pkin(8) + t178 * t167 + pkin(2);
	t160 = t167 * pkin(8) - t178 * t165;
	t171 = sin(qJ(2));
	t175 = cos(qJ(2));
	t179 = t159 * t171 - t160 * t175;
	t177 = t170 * pkin(4) - pkin(9) * t174 + pkin(7) + qJ(3);
	t163 = cos(t164);
	t161 = t166 * t162 * t174 + t168 * t170;
	t158 = t159 * t175 + t160 * t171 + pkin(1);
	t157 = t166 * t177 - t179 * t168;
	t1 = [(t168 * t191 + t180) * t163 + (-t168 * t183 + t190) * t162 + t172 * t184, (t168 * t188 - t181) * t163 + (t168 * t182 + t187) * t162 - t172 * t185, -(-t176 * t163 + t172 * t194) * t170 - t166 * t189, t157 * t172 + t158 * t176 + 0; (-t168 * t190 + t183) * t163 + (t168 * t180 + t191) * t162 - t176 * t184, (-t168 * t187 - t182) * t163 + (-t168 * t181 + t188) * t162 + t176 * t185, (t172 * t163 + t176 * t194) * t170 + t166 * t186, -t157 * t176 + t158 * t172 + 0; -t166 * t163 * t169 + t161 * t173, -t161 * t169 - t163 * t192, t162 * t193 - t168 * t174, t179 * t166 + t168 * t177 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end