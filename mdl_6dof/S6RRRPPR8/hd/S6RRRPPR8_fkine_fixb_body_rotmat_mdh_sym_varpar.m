% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRPPR8 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:24
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRRPPR8_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPR8_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:24:51
	% EndTime: 2020-11-04 22:24:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:24:51
	% EndTime: 2020-11-04 22:24:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t97 = cos(qJ(1));
	t96 = sin(qJ(1));
	t1 = [t97, -t96, 0, 0; t96, t97, 0, 0; 0, 0, 1, pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:24:51
	% EndTime: 2020-11-04 22:24:51
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t101 = sin(qJ(1));
	t98 = sin(pkin(6));
	t109 = t101 * t98;
	t103 = cos(qJ(1));
	t108 = t103 * t98;
	t100 = sin(qJ(2));
	t107 = t101 * t100;
	t102 = cos(qJ(2));
	t106 = t101 * t102;
	t105 = t103 * t100;
	t104 = t103 * t102;
	t99 = cos(pkin(6));
	t1 = [-t99 * t107 + t104, -t99 * t106 - t105, t109, t103 * pkin(1) + pkin(8) * t109 + 0; t99 * t105 + t106, t99 * t104 - t107, -t108, t101 * pkin(1) - pkin(8) * t108 + 0; t98 * t100, t98 * t102, t99, t99 * pkin(8) + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:24:51
	% EndTime: 2020-11-04 22:24:51
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->23), mult. (60->37), div. (0->0), fcn. (81->8), ass. (0->20)
	t113 = sin(pkin(6));
	t118 = cos(qJ(3));
	t128 = t113 * t118;
	t115 = sin(qJ(3));
	t127 = t115 * t113;
	t116 = sin(qJ(2));
	t126 = t116 * t118;
	t117 = sin(qJ(1));
	t125 = t117 * t116;
	t119 = cos(qJ(2));
	t124 = t117 * t119;
	t120 = cos(qJ(1));
	t123 = t120 * t116;
	t122 = t120 * t119;
	t121 = pkin(2) * t116 - pkin(9) * t119;
	t114 = cos(pkin(6));
	t112 = t119 * pkin(2) + t116 * pkin(9) + pkin(1);
	t111 = -t114 * t123 - t124;
	t110 = t113 * pkin(8) - t121 * t114;
	t1 = [(-t114 * t126 + t127) * t117 + t118 * t122, (t114 * t125 - t122) * t115 + t117 * t128, t114 * t124 + t123, t110 * t117 + t112 * t120 + 0; -t111 * t118 - t120 * t127, t111 * t115 - t120 * t128, -t114 * t122 + t125, -t110 * t120 + t112 * t117 + 0; t113 * t126 + t114 * t115, t114 * t118 - t116 * t127, -t113 * t119, t114 * pkin(8) + t121 * t113 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:24:51
	% EndTime: 2020-11-04 22:24:51
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (45->27), mult. (78->41), div. (0->0), fcn. (99->8), ass. (0->22)
	t133 = sin(pkin(6));
	t138 = cos(qJ(3));
	t149 = t133 * t138;
	t135 = sin(qJ(3));
	t148 = t135 * t133;
	t136 = sin(qJ(2));
	t147 = t136 * t138;
	t137 = sin(qJ(1));
	t146 = t137 * t136;
	t139 = cos(qJ(2));
	t145 = t137 * t139;
	t140 = cos(qJ(1));
	t144 = t140 * t136;
	t143 = t140 * t139;
	t132 = pkin(3) * t138 + qJ(4) * t135 + pkin(2);
	t142 = pkin(9) * t139 - t132 * t136;
	t141 = pkin(3) * t135 - qJ(4) * t138 + pkin(8);
	t134 = cos(pkin(6));
	t131 = t134 * t144 + t145;
	t130 = pkin(9) * t136 + t132 * t139 + pkin(1);
	t129 = t133 * t141 + t134 * t142;
	t1 = [(-t134 * t147 + t148) * t137 + t138 * t143, t134 * t145 + t144, (-t134 * t146 + t143) * t135 - t137 * t149, t129 * t137 + t130 * t140 + 0; t131 * t138 - t140 * t148, -t134 * t143 + t146, t131 * t135 + t140 * t149, -t129 * t140 + t130 * t137 + 0; t133 * t147 + t134 * t135, -t133 * t139, -t134 * t138 + t136 * t148, -t133 * t142 + t134 * t141 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:24:51
	% EndTime: 2020-11-04 22:24:52
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (57->28), mult. (78->41), div. (0->0), fcn. (99->8), ass. (0->24)
	t154 = sin(pkin(6));
	t160 = cos(qJ(3));
	t172 = t154 * t160;
	t157 = sin(qJ(3));
	t171 = t157 * t154;
	t158 = sin(qJ(2));
	t170 = t158 * t160;
	t159 = sin(qJ(1));
	t169 = t159 * t158;
	t161 = cos(qJ(2));
	t168 = t159 * t161;
	t162 = cos(qJ(1));
	t167 = t162 * t158;
	t166 = t162 * t161;
	t163 = pkin(3) + pkin(4);
	t153 = qJ(4) * t157 + t163 * t160 + pkin(2);
	t156 = qJ(5) - pkin(9);
	t165 = t153 * t158 + t156 * t161;
	t164 = qJ(4) * t160 - t163 * t157 - pkin(8);
	t155 = cos(pkin(6));
	t152 = -t155 * t167 - t168;
	t151 = t153 * t161 - t156 * t158 + pkin(1);
	t150 = t154 * t164 + t165 * t155;
	t1 = [(-t155 * t169 + t166) * t157 - t159 * t172, (t155 * t170 - t171) * t159 - t160 * t166, -t155 * t168 - t167, -t150 * t159 + t151 * t162 + 0; -t152 * t157 + t162 * t172, t152 * t160 + t162 * t171, t155 * t166 - t169, t150 * t162 + t151 * t159 + 0; -t155 * t160 + t158 * t171, -t154 * t170 - t155 * t157, t154 * t161, t165 * t154 - t164 * t155 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:24:52
	% EndTime: 2020-11-04 22:24:52
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (86->40), mult. (116->63), div. (0->0), fcn. (150->10), ass. (0->33)
	t182 = pkin(3) + pkin(4) + pkin(10);
	t188 = sin(qJ(3));
	t179 = t182 * t188 + pkin(8);
	t186 = qJ(4) + pkin(5);
	t180 = t186 * t188 + pkin(2);
	t183 = sin(pkin(6));
	t184 = cos(pkin(6));
	t189 = sin(qJ(2));
	t192 = cos(qJ(3));
	t185 = qJ(5) - pkin(9);
	t193 = cos(qJ(2));
	t206 = t185 * t193;
	t208 = (t180 * t189 + t206) * t184 + (t184 * t182 * t189 + t183 * t186) * t192 - t179 * t183;
	t207 = t184 * t192;
	t187 = sin(qJ(6));
	t205 = t187 * t193;
	t204 = t188 * t183;
	t203 = t188 * t189;
	t190 = sin(qJ(1));
	t202 = t190 * t193;
	t191 = cos(qJ(6));
	t201 = t191 * t193;
	t200 = t192 * t183;
	t194 = cos(qJ(1));
	t199 = t193 * t194;
	t198 = t194 * t189;
	t175 = t184 * t203 + t200;
	t195 = -t175 * t187 + t184 * t201;
	t178 = t182 * t192 + t180;
	t177 = t188 * t205 + t191 * t189;
	t176 = t183 * t203 - t207;
	t173 = t178 * t193 - t185 * t189 + pkin(1);
	t1 = [(-t175 * t190 + t188 * t199) * t191 - t187 * (t184 * t202 + t198), -t194 * t177 - t195 * t190, (-t189 * t207 + t204) * t190 + t192 * t199, t173 * t194 - t208 * t190 + 0; (t175 * t191 + t184 * t205) * t194 + t190 * (-t187 * t189 + t188 * t201), -t190 * t177 + t195 * t194, (t184 * t198 + t202) * t192 - t194 * t204, t173 * t190 + t208 * t194 + 0; t176 * t191 + t183 * t205, -t176 * t187 + t183 * t201, t184 * t188 + t189 * t200, pkin(7) + 0 + (t178 * t189 + t206) * t183 + (-t186 * t192 + t179) * t184; 0, 0, 0, 1;];
	Tc_mdh = t1;
end