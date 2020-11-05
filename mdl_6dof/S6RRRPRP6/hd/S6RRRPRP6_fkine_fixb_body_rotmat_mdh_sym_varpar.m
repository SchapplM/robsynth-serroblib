% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRPRP6 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:27
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRRPRP6_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP6_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRP6_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP6_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:27:48
	% EndTime: 2020-11-04 22:27:48
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:27:48
	% EndTime: 2020-11-04 22:27:48
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t102 = cos(qJ(1));
	t101 = sin(qJ(1));
	t1 = [t102, -t101, 0, 0; t101, t102, 0, 0; 0, 0, 1, pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:27:48
	% EndTime: 2020-11-04 22:27:48
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t103 = sin(pkin(6));
	t106 = sin(qJ(1));
	t114 = t106 * t103;
	t105 = sin(qJ(2));
	t113 = t106 * t105;
	t107 = cos(qJ(2));
	t112 = t106 * t107;
	t108 = cos(qJ(1));
	t111 = t108 * t103;
	t110 = t108 * t105;
	t109 = t108 * t107;
	t104 = cos(pkin(6));
	t1 = [-t104 * t113 + t109, -t104 * t112 - t110, t114, t108 * pkin(1) + pkin(8) * t114 + 0; t104 * t110 + t112, t104 * t109 - t113, -t111, t106 * pkin(1) - pkin(8) * t111 + 0; t103 * t105, t103 * t107, t104, t104 * pkin(8) + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:27:48
	% EndTime: 2020-11-04 22:27:49
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (29->23), mult. (60->37), div. (0->0), fcn. (81->8), ass. (0->20)
	t118 = sin(pkin(6));
	t123 = cos(qJ(3));
	t133 = t118 * t123;
	t120 = sin(qJ(3));
	t132 = t120 * t118;
	t121 = sin(qJ(2));
	t131 = t121 * t123;
	t122 = sin(qJ(1));
	t130 = t122 * t121;
	t124 = cos(qJ(2));
	t129 = t122 * t124;
	t125 = cos(qJ(1));
	t128 = t125 * t121;
	t127 = t125 * t124;
	t126 = pkin(2) * t121 - pkin(9) * t124;
	t119 = cos(pkin(6));
	t117 = t124 * pkin(2) + t121 * pkin(9) + pkin(1);
	t116 = -t119 * t128 - t129;
	t115 = t118 * pkin(8) - t126 * t119;
	t1 = [(-t119 * t131 + t132) * t122 + t123 * t127, (t119 * t130 - t127) * t120 + t122 * t133, t119 * t129 + t128, t115 * t122 + t117 * t125 + 0; -t116 * t123 - t125 * t132, t116 * t120 - t125 * t133, -t119 * t127 + t130, -t115 * t125 + t117 * t122 + 0; t118 * t131 + t119 * t120, t119 * t123 - t121 * t132, -t118 * t124, t119 * pkin(8) + t126 * t118 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:27:49
	% EndTime: 2020-11-04 22:27:49
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (54->26), mult. (69->38), div. (0->0), fcn. (90->10), ass. (0->25)
	t142 = sin(pkin(6));
	t146 = sin(qJ(2));
	t158 = t142 * t146;
	t147 = sin(qJ(1));
	t157 = t142 * t147;
	t149 = cos(qJ(1));
	t156 = t142 * t149;
	t155 = t147 * t146;
	t148 = cos(qJ(2));
	t154 = t147 * t148;
	t153 = t149 * t146;
	t152 = t149 * t148;
	t151 = sin(qJ(3)) * pkin(3) + pkin(8);
	t138 = cos(qJ(3)) * pkin(3) + pkin(2);
	t144 = qJ(4) + pkin(9);
	t150 = t138 * t146 - t144 * t148;
	t143 = cos(pkin(6));
	t141 = qJ(3) + pkin(11);
	t140 = cos(t141);
	t139 = sin(t141);
	t137 = t143 * t153 + t154;
	t136 = t143 * t155 - t152;
	t135 = t138 * t148 + t144 * t146 + pkin(1);
	t134 = t142 * t151 - t150 * t143;
	t1 = [-t136 * t140 + t139 * t157, t136 * t139 + t140 * t157, t143 * t154 + t153, t134 * t147 + t135 * t149 + 0; t137 * t140 - t139 * t156, -t137 * t139 - t140 * t156, -t143 * t152 + t155, -t134 * t149 + t135 * t147 + 0; t143 * t139 + t140 * t158, -t139 * t158 + t143 * t140, -t142 * t148, t150 * t142 + t151 * t143 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:27:49
	% EndTime: 2020-11-04 22:27:49
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (104->39), mult. (152->61), div. (0->0), fcn. (186->14), ass. (0->40)
	t176 = sin(pkin(11));
	t178 = cos(pkin(11));
	t170 = t178 * pkin(4) + t176 * pkin(10) + pkin(3);
	t182 = sin(qJ(3));
	t167 = t170 * t182 + pkin(8);
	t171 = -t176 * pkin(4) + t178 * pkin(10);
	t168 = t171 * t182 + pkin(2);
	t177 = sin(pkin(6));
	t179 = cos(pkin(6));
	t183 = sin(qJ(2));
	t186 = cos(qJ(3));
	t180 = qJ(4) + pkin(9);
	t187 = cos(qJ(2));
	t197 = t180 * t187;
	t202 = (t168 * t183 - t197) * t179 + (t179 * t170 * t183 + t177 * t171) * t186 - t177 * t167;
	t201 = t177 * t183;
	t184 = sin(qJ(1));
	t200 = t177 * t184;
	t199 = t177 * t187;
	t188 = cos(qJ(1));
	t198 = t177 * t188;
	t196 = t184 * t183;
	t195 = t184 * t187;
	t194 = t188 * t183;
	t193 = t188 * t187;
	t163 = t179 * t196 - t193;
	t175 = qJ(3) + pkin(11);
	t173 = sin(t175);
	t174 = cos(t175);
	t190 = -t163 * t174 + t173 * t200;
	t165 = t179 * t194 + t195;
	t189 = -t165 * t174 + t173 * t198;
	t185 = cos(qJ(5));
	t181 = sin(qJ(5));
	t166 = t179 * t195 + t194;
	t164 = t179 * t193 - t196;
	t162 = t179 * t173 + t174 * t201;
	t161 = t170 * t186 + t168;
	t159 = t161 * t187 + t180 * t183 + pkin(1);
	t1 = [t166 * t181 + t190 * t185, t166 * t185 - t190 * t181, -t163 * t173 - t174 * t200, t159 * t188 - t202 * t184 + 0; -t181 * t164 - t189 * t185, -t185 * t164 + t189 * t181, t165 * t173 + t174 * t198, t159 * t184 + t202 * t188 + 0; t162 * t185 - t181 * t199, -t162 * t181 - t185 * t199, t173 * t201 - t179 * t174, pkin(7) + 0 + (t161 * t183 - t197) * t177 + (-t171 * t186 + t167) * t179; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:27:49
	% EndTime: 2020-11-04 22:27:49
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (111->44), mult. (155->61), div. (0->0), fcn. (201->12), ass. (0->40)
	t221 = cos(pkin(6));
	t229 = cos(qJ(2));
	t230 = cos(qJ(1));
	t232 = t230 * t229;
	t226 = sin(qJ(2));
	t227 = sin(qJ(1));
	t235 = t227 * t226;
	t212 = -t221 * t232 + t235;
	t224 = sin(qJ(5));
	t242 = t212 * t224;
	t233 = t230 * t226;
	t234 = t227 * t229;
	t214 = t221 * t234 + t233;
	t241 = t214 * t224;
	t216 = cos(qJ(3)) * pkin(3) + pkin(2);
	t240 = t216 * t226;
	t220 = sin(pkin(6));
	t239 = t220 * t226;
	t238 = t220 * t227;
	t237 = t220 * t229;
	t236 = t220 * t230;
	t231 = sin(qJ(3)) * pkin(3) + pkin(8);
	t228 = cos(qJ(5));
	t223 = qJ(4) + pkin(9);
	t222 = -qJ(6) - pkin(10);
	t219 = qJ(3) + pkin(11);
	t218 = cos(t219);
	t217 = sin(t219);
	t215 = t228 * pkin(5) + pkin(4);
	t213 = t221 * t233 + t234;
	t211 = t221 * t235 - t232;
	t210 = t216 * t229 + t223 * t226 + pkin(1);
	t209 = t221 * t217 + t218 * t239;
	t208 = t217 * t239 - t221 * t218;
	t207 = t220 * t231 + (t223 * t229 - t240) * t221;
	t206 = -t211 * t218 + t217 * t238;
	t205 = t213 * t218 - t217 * t236;
	t204 = t213 * t217 + t218 * t236;
	t203 = t211 * t217 + t218 * t238;
	t1 = [t206 * t228 + t241, -t206 * t224 + t214 * t228, -t203, pkin(5) * t241 + t203 * t222 + t206 * t215 + t207 * t227 + t210 * t230 + 0; t205 * t228 + t242, -t205 * t224 + t212 * t228, t204, pkin(5) * t242 - t204 * t222 + t205 * t215 - t207 * t230 + t210 * t227 + 0; t209 * t228 - t224 * t237, -t209 * t224 - t228 * t237, t208, -t208 * t222 + t209 * t215 + pkin(7) + 0 + t231 * t221 + (t240 + (-pkin(5) * t224 - t223) * t229) * t220; 0, 0, 0, 1;];
	Tc_mdh = t1;
end