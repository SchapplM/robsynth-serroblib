% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRRP10 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:45
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRRRRP10_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP10_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRP10_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP10_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:45:19
	% EndTime: 2020-11-04 22:45:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:45:19
	% EndTime: 2020-11-04 22:45:19
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t96 = cos(qJ(1));
	t95 = sin(qJ(1));
	t1 = [t96, -t95, 0, 0; t95, t96, 0, 0; 0, 0, 1, pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:45:19
	% EndTime: 2020-11-04 22:45:19
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t100 = sin(qJ(1));
	t97 = sin(pkin(6));
	t108 = t100 * t97;
	t99 = sin(qJ(2));
	t107 = t100 * t99;
	t102 = cos(qJ(1));
	t106 = t102 * t97;
	t105 = t102 * t99;
	t101 = cos(qJ(2));
	t104 = t100 * t101;
	t103 = t102 * t101;
	t98 = cos(pkin(6));
	t1 = [-t98 * t107 + t103, -t98 * t104 - t105, t108, t102 * pkin(1) + pkin(8) * t108 + 0; t98 * t105 + t104, t98 * t103 - t107, -t106, t100 * pkin(1) - pkin(8) * t106 + 0; t97 * t99, t97 * t101, t98, t98 * pkin(8) + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:45:19
	% EndTime: 2020-11-04 22:45:19
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->23), mult. (60->37), div. (0->0), fcn. (81->8), ass. (0->20)
	t112 = sin(pkin(6));
	t117 = cos(qJ(3));
	t127 = t112 * t117;
	t114 = sin(qJ(3));
	t126 = t114 * t112;
	t115 = sin(qJ(2));
	t125 = t115 * t117;
	t116 = sin(qJ(1));
	t124 = t116 * t115;
	t118 = cos(qJ(2));
	t123 = t116 * t118;
	t119 = cos(qJ(1));
	t122 = t119 * t115;
	t121 = t119 * t118;
	t120 = pkin(2) * t115 - pkin(9) * t118;
	t113 = cos(pkin(6));
	t111 = t118 * pkin(2) + t115 * pkin(9) + pkin(1);
	t110 = t113 * t122 + t123;
	t109 = t112 * pkin(8) - t120 * t113;
	t1 = [(-t113 * t125 + t126) * t116 + t117 * t121, (t113 * t124 - t121) * t114 + t116 * t127, t113 * t123 + t122, t109 * t116 + t111 * t119 + 0; t110 * t117 - t119 * t126, -t110 * t114 - t119 * t127, -t113 * t121 + t124, -t109 * t119 + t111 * t116 + 0; t112 * t125 + t113 * t114, t113 * t117 - t115 * t126, -t112 * t118, t113 * pkin(8) + t120 * t112 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:45:19
	% EndTime: 2020-11-04 22:45:19
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (57->34), mult. (112->60), div. (0->0), fcn. (146->10), ass. (0->29)
	t135 = sin(pkin(6));
	t142 = cos(qJ(3));
	t156 = t135 * t142;
	t137 = sin(qJ(4));
	t143 = cos(qJ(2));
	t155 = t137 * t143;
	t138 = sin(qJ(3));
	t154 = t138 * t135;
	t139 = sin(qJ(2));
	t153 = t139 * t142;
	t140 = sin(qJ(1));
	t152 = t140 * t143;
	t141 = cos(qJ(4));
	t151 = t141 * t143;
	t150 = t142 * t143;
	t144 = cos(qJ(1));
	t149 = t144 * t139;
	t148 = t144 * t143;
	t133 = t142 * pkin(3) + t138 * pkin(10) + pkin(2);
	t147 = t143 * pkin(9) - t133 * t139;
	t146 = t138 * pkin(3) - t142 * pkin(10) + pkin(8);
	t136 = cos(pkin(6));
	t131 = t136 * t153 - t154;
	t145 = t131 * t137 + t136 * t151;
	t132 = t137 * t150 - t141 * t139;
	t130 = t135 * t153 + t136 * t138;
	t129 = t139 * pkin(9) + t133 * t143 + pkin(1);
	t128 = t135 * t146 + t147 * t136;
	t1 = [(-t131 * t140 + t142 * t148) * t141 + (t136 * t152 + t149) * t137, -t144 * t132 + t145 * t140, (-t140 * t136 * t139 + t148) * t138 - t140 * t156, t128 * t140 + t129 * t144 + 0; (t131 * t141 - t136 * t155) * t144 + t140 * (t137 * t139 + t141 * t150), -t140 * t132 - t145 * t144, (t136 * t149 + t152) * t138 + t144 * t156, -t128 * t144 + t129 * t140 + 0; t130 * t141 - t135 * t155, -t130 * t137 - t135 * t151, -t136 * t142 + t139 * t154, -t147 * t135 + t146 * t136 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:45:19
	% EndTime: 2020-11-04 22:45:19
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (90->36), mult. (125->55), div. (0->0), fcn. (159->12), ass. (0->34)
	t170 = sin(pkin(6));
	t175 = cos(qJ(3));
	t190 = t170 * t175;
	t176 = cos(qJ(2));
	t189 = t170 * t176;
	t172 = sin(qJ(3));
	t188 = t172 * t170;
	t173 = sin(qJ(2));
	t187 = t173 * t175;
	t174 = sin(qJ(1));
	t186 = t174 * t173;
	t185 = t174 * t176;
	t177 = cos(qJ(1));
	t184 = t177 * t173;
	t183 = t177 * t176;
	t166 = cos(qJ(4)) * pkin(4) + pkin(3);
	t178 = pkin(11) + pkin(10);
	t159 = t166 * t175 + t178 * t172 + pkin(2);
	t165 = sin(qJ(4)) * pkin(4) + pkin(9);
	t182 = t159 * t173 - t165 * t176;
	t181 = t166 * t172 - t178 * t175 + pkin(8);
	t171 = cos(pkin(6));
	t161 = t171 * t187 - t188;
	t180 = -t161 * t174 + t175 * t183;
	t179 = t161 * t177 + t175 * t185;
	t169 = qJ(4) + qJ(5);
	t168 = cos(t169);
	t167 = sin(t169);
	t163 = t171 * t185 + t184;
	t162 = t171 * t183 - t186;
	t160 = t170 * t187 + t171 * t172;
	t158 = t159 * t176 + t165 * t173 + pkin(1);
	t157 = t181 * t170 - t182 * t171;
	t1 = [t163 * t167 + t180 * t168, t163 * t168 - t180 * t167, (-t171 * t186 + t183) * t172 - t174 * t190, t157 * t174 + t158 * t177 + 0; -t167 * t162 + t179 * t168, -t168 * t162 - t179 * t167, (t171 * t184 + t185) * t172 + t177 * t190, -t157 * t177 + t158 * t174 + 0; t160 * t168 - t167 * t189, -t160 * t167 - t168 * t189, -t171 * t175 + t173 * t188, t182 * t170 + t181 * t171 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:45:19
	% EndTime: 2020-11-04 22:45:19
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (128->44), mult. (187->61), div. (0->0), fcn. (241->12), ass. (0->40)
	t210 = sin(pkin(6));
	t215 = cos(qJ(3));
	t230 = t210 * t215;
	t216 = cos(qJ(2));
	t229 = t210 * t216;
	t212 = sin(qJ(3));
	t228 = t212 * t210;
	t213 = sin(qJ(2));
	t227 = t213 * t215;
	t214 = sin(qJ(1));
	t226 = t214 * t213;
	t225 = t214 * t216;
	t217 = cos(qJ(1));
	t224 = t217 * t213;
	t223 = t217 * t216;
	t206 = cos(qJ(4)) * pkin(4) + pkin(3);
	t218 = pkin(11) + pkin(10);
	t199 = t206 * t215 + t218 * t212 + pkin(2);
	t205 = sin(qJ(4)) * pkin(4) + pkin(9);
	t222 = t199 * t213 - t205 * t216;
	t221 = t206 * t212 - t218 * t215 + pkin(8);
	t211 = cos(pkin(6));
	t201 = t211 * t227 - t228;
	t220 = -t201 * t214 + t215 * t223;
	t219 = t201 * t217 + t215 * t225;
	t209 = qJ(4) + qJ(5);
	t208 = cos(t209);
	t207 = sin(t209);
	t203 = t211 * t225 + t224;
	t202 = t211 * t223 - t226;
	t200 = t210 * t227 + t211 * t212;
	t198 = t199 * t216 + t205 * t213 + pkin(1);
	t197 = t200 * t208 - t207 * t229;
	t196 = t200 * t207 + t208 * t229;
	t195 = -t208 * t202 - t219 * t207;
	t194 = t203 * t207 + t220 * t208;
	t193 = -t207 * t202 + t219 * t208;
	t192 = t203 * t208 - t220 * t207;
	t191 = t221 * t210 - t222 * t211;
	t1 = [t194, (-t211 * t226 + t223) * t212 - t214 * t230, -t192, t194 * pkin(5) - t192 * qJ(6) + t191 * t214 + t198 * t217 + 0; t193, (t211 * t224 + t225) * t212 + t217 * t230, -t195, t193 * pkin(5) - t195 * qJ(6) - t191 * t217 + t198 * t214 + 0; t197, -t211 * t215 + t213 * t228, t196, t197 * pkin(5) + t196 * qJ(6) + t222 * t210 + t221 * t211 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end