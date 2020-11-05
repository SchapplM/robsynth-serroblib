% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRPPRR3 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:58
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PRPPRR3_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPPRR3_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:58:11
	% EndTime: 2020-11-04 20:58:11
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:58:11
	% EndTime: 2020-11-04 20:58:11
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t97 = cos(pkin(10));
	t96 = sin(pkin(10));
	t1 = [t97, -t96, 0, 0; t96, t97, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:58:11
	% EndTime: 2020-11-04 20:58:11
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t98 = sin(pkin(10));
	t99 = sin(pkin(6));
	t107 = t98 * t99;
	t100 = cos(pkin(10));
	t106 = t100 * t99;
	t101 = cos(pkin(6));
	t102 = sin(qJ(2));
	t105 = t101 * t102;
	t103 = cos(qJ(2));
	t104 = t101 * t103;
	t1 = [t100 * t103 - t98 * t105, -t100 * t102 - t98 * t104, t107, t100 * pkin(1) + pkin(7) * t107 + 0; t100 * t105 + t98 * t103, t100 * t104 - t98 * t102, -t106, t98 * pkin(1) - pkin(7) * t106 + 0; t99 * t102, t99 * t103, t101, t101 * pkin(7) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:58:11
	% EndTime: 2020-11-04 20:58:11
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (22->22), mult. (43->36), div. (0->0), fcn. (56->6), ass. (0->13)
	t108 = sin(pkin(10));
	t109 = sin(pkin(6));
	t119 = t108 * t109;
	t111 = cos(pkin(6));
	t118 = t108 * t111;
	t110 = cos(pkin(10));
	t117 = t110 * t109;
	t116 = t110 * t111;
	t112 = sin(qJ(2));
	t115 = t111 * t112;
	t113 = cos(qJ(2));
	t114 = t111 * t113;
	t1 = [-t108 * t115 + t110 * t113, t119, t108 * t114 + t110 * t112, (t110 * pkin(2) + qJ(3) * t118) * t113 + (-pkin(2) * t118 + t110 * qJ(3)) * t112 + pkin(7) * t119 + t110 * pkin(1) + 0; t108 * t113 + t110 * t115, -t117, t108 * t112 - t110 * t114, (t108 * pkin(2) - qJ(3) * t116) * t113 + (pkin(2) * t116 + t108 * qJ(3)) * t112 - pkin(7) * t117 + t108 * pkin(1) + 0; t109 * t112, t111, -t109 * t113, t111 * pkin(7) + qJ(1) + 0 + (pkin(2) * t112 - qJ(3) * t113) * t109; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:58:11
	% EndTime: 2020-11-04 20:58:11
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (41->31), mult. (67->47), div. (0->0), fcn. (90->8), ass. (0->21)
	t125 = sin(pkin(10));
	t126 = sin(pkin(6));
	t139 = t125 * t126;
	t129 = cos(pkin(6));
	t138 = t125 * t129;
	t127 = cos(pkin(11));
	t137 = t127 * t125;
	t128 = cos(pkin(10));
	t136 = t128 * t126;
	t135 = t128 * t129;
	t133 = pkin(2) + pkin(3);
	t134 = t129 * t133;
	t132 = cos(qJ(2));
	t131 = sin(qJ(2));
	t130 = qJ(4) - pkin(7);
	t124 = sin(pkin(11));
	t123 = t124 * t125 + t127 * t135;
	t122 = t124 * t128 - t129 * t137;
	t121 = t124 * t135 - t137;
	t120 = t124 * t138 + t128 * t127;
	t1 = [t120 * t132 + t131 * t122, t131 * t120 - t132 * t122, -t139, (qJ(3) * t138 + t128 * t133) * t132 + (t128 * qJ(3) - t125 * t134) * t131 - t130 * t139 + t128 * pkin(1) + 0; -t132 * t121 + t131 * t123, -t121 * t131 - t132 * t123, t136, (-qJ(3) * t135 + t125 * t133) * t132 + (t125 * qJ(3) + t128 * t134) * t131 + t130 * t136 + t125 * pkin(1) + 0; t126 * (-t124 * t132 + t127 * t131), -t126 * (t124 * t131 + t127 * t132), -t129, -t130 * t129 + qJ(1) + 0 + (-qJ(3) * t132 + t131 * t133) * t126; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:58:11
	% EndTime: 2020-11-04 20:58:11
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (71->38), mult. (120->62), div. (0->0), fcn. (156->10), ass. (0->27)
	t148 = sin(pkin(10));
	t152 = cos(pkin(6));
	t165 = t148 * t152;
	t149 = sin(pkin(6));
	t153 = qJ(4) - pkin(7);
	t164 = t149 * t153;
	t154 = sin(qJ(5));
	t163 = t149 * t154;
	t156 = cos(qJ(5));
	t162 = t149 * t156;
	t150 = cos(pkin(11));
	t161 = t150 * t148;
	t151 = cos(pkin(10));
	t160 = t151 * t152;
	t147 = sin(pkin(11));
	t140 = t147 * t165 + t151 * t150;
	t142 = -t147 * t151 + t152 * t161;
	t155 = sin(qJ(2));
	t157 = cos(qJ(2));
	t159 = t140 * t157 - t155 * t142;
	t141 = t147 * t160 - t161;
	t143 = t147 * t148 + t150 * t160;
	t158 = t157 * t141 - t155 * t143;
	t146 = t147 * pkin(4) - t150 * pkin(8) + qJ(3);
	t145 = t150 * pkin(4) + t147 * pkin(8) + pkin(2) + pkin(3);
	t144 = t147 * t157 - t150 * t155;
	t1 = [-t148 * t163 + t159 * t156, -t148 * t162 - t159 * t154, -t140 * t155 - t157 * t142, (t151 * t145 + t146 * t165) * t157 + (-t145 * t165 + t151 * t146) * t155 - t148 * t164 + t151 * pkin(1) + 0; t151 * t163 - t158 * t156, t151 * t162 + t158 * t154, t155 * t141 + t157 * t143, (t148 * t145 - t146 * t160) * t157 + (t145 * t160 + t148 * t146) * t155 + t151 * t164 + t148 * pkin(1) + 0; -t144 * t162 - t152 * t154, t144 * t163 - t152 * t156, t149 * (t147 * t155 + t150 * t157), -t153 * t152 + qJ(1) + 0 + (t145 * t155 - t146 * t157) * t149; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:58:11
	% EndTime: 2020-11-04 20:58:12
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (120->53), mult. (238->80), div. (0->0), fcn. (317->12), ass. (0->38)
	t182 = sin(pkin(11));
	t184 = sin(pkin(6));
	t185 = cos(pkin(11));
	t191 = sin(qJ(2));
	t194 = cos(qJ(2));
	t203 = (t182 * t191 + t185 * t194) * t184;
	t183 = sin(pkin(10));
	t187 = cos(pkin(6));
	t202 = t183 * t187;
	t188 = qJ(4) - pkin(7);
	t201 = t184 * t188;
	t190 = sin(qJ(5));
	t200 = t184 * t190;
	t193 = cos(qJ(5));
	t199 = t184 * t193;
	t198 = t185 * t183;
	t186 = cos(pkin(10));
	t197 = t186 * t187;
	t174 = t182 * t202 + t186 * t185;
	t176 = -t182 * t186 + t187 * t198;
	t196 = t174 * t194 - t191 * t176;
	t175 = t182 * t197 - t198;
	t177 = t182 * t183 + t185 * t197;
	t195 = t194 * t175 - t191 * t177;
	t192 = cos(qJ(6));
	t189 = sin(qJ(6));
	t181 = t182 * pkin(4) - t185 * pkin(8) + qJ(3);
	t180 = t185 * pkin(4) + t182 * pkin(8) + pkin(2) + pkin(3);
	t178 = t182 * t194 - t185 * t191;
	t173 = -t178 * t199 - t187 * t190;
	t172 = t178 * t200 - t187 * t193;
	t171 = -t174 * t191 - t194 * t176;
	t170 = t191 * t175 + t194 * t177;
	t169 = t186 * t200 - t195 * t193;
	t168 = -t183 * t199 - t196 * t190;
	t167 = -t183 * t200 + t196 * t193;
	t166 = t186 * t199 + t195 * t190;
	t1 = [t167 * t192 + t171 * t189, -t167 * t189 + t171 * t192, -t168, t167 * pkin(5) - t168 * pkin(9) + (t186 * t180 + t181 * t202) * t194 + (-t180 * t202 + t186 * t181) * t191 - t183 * t201 + t186 * pkin(1) + 0; t169 * t192 + t170 * t189, -t169 * t189 + t170 * t192, -t166, t169 * pkin(5) - t166 * pkin(9) + (t183 * t180 - t181 * t197) * t194 + (t180 * t197 + t183 * t181) * t191 + t186 * t201 + t183 * pkin(1) + 0; t173 * t192 + t189 * t203, -t173 * t189 + t192 * t203, -t172, t173 * pkin(5) - t172 * pkin(9) - t188 * t187 + qJ(1) + 0 + (t180 * t191 - t181 * t194) * t184; 0, 0, 0, 1;];
	Tc_mdh = t1;
end