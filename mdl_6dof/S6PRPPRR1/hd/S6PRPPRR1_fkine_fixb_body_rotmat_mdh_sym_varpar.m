% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRPPRR1 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:57
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PRPPRR1_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPPRR1_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:57:17
	% EndTime: 2020-11-04 20:57:17
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:57:17
	% EndTime: 2020-11-04 20:57:17
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t100 = cos(pkin(10));
	t99 = sin(pkin(10));
	t1 = [t100, -t99, 0, 0; t99, t100, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:57:17
	% EndTime: 2020-11-04 20:57:17
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t101 = sin(pkin(10));
	t102 = sin(pkin(6));
	t110 = t101 * t102;
	t103 = cos(pkin(10));
	t109 = t103 * t102;
	t104 = cos(pkin(6));
	t105 = sin(qJ(2));
	t108 = t104 * t105;
	t106 = cos(qJ(2));
	t107 = t104 * t106;
	t1 = [-t101 * t108 + t103 * t106, -t101 * t107 - t103 * t105, t110, t103 * pkin(1) + pkin(7) * t110 + 0; t101 * t106 + t103 * t108, -t101 * t105 + t103 * t107, -t109, t101 * pkin(1) - pkin(7) * t109 + 0; t102 * t105, t102 * t106, t104, t104 * pkin(7) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:57:17
	% EndTime: 2020-11-04 20:57:17
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (53->22), mult. (37->26), div. (0->0), fcn. (44->12), ass. (0->20)
	t130 = pkin(2) * sin(qJ(2));
	t123 = qJ(2) + pkin(11);
	t128 = pkin(7) + qJ(3);
	t127 = cos(pkin(6));
	t126 = cos(pkin(10));
	t125 = sin(pkin(6));
	t124 = sin(pkin(10));
	t122 = cos(t123);
	t121 = sin(t123);
	t120 = pkin(6) - t123;
	t119 = pkin(6) + t123;
	t118 = cos(qJ(2)) * pkin(2) + pkin(1);
	t117 = cos(t119);
	t116 = sin(t120);
	t115 = cos(t120) / 0.2e1;
	t114 = sin(t119) / 0.2e1;
	t113 = -t125 * t128 + t127 * t130;
	t112 = t117 / 0.2e1 + t115;
	t111 = t114 - t116 / 0.2e1;
	t1 = [-t124 * t111 + t126 * t122, -t124 * t112 - t126 * t121, t124 * t125, -t124 * t113 + t126 * t118 + 0; t126 * t111 + t124 * t122, t126 * t112 - t124 * t121, -t126 * t125, t126 * t113 + t124 * t118 + 0; t115 - t117 / 0.2e1, t116 / 0.2e1 + t114, t127, t125 * t130 + t127 * t128 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:57:17
	% EndTime: 2020-11-04 20:57:17
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (73->35), mult. (85->54), div. (0->0), fcn. (105->16), ass. (0->28)
	t141 = sin(pkin(12));
	t144 = sin(pkin(6));
	t157 = t141 * t144;
	t143 = sin(pkin(10));
	t148 = cos(pkin(6));
	t156 = t143 * t148;
	t145 = cos(pkin(12));
	t155 = t144 * t145;
	t147 = cos(pkin(10));
	t154 = t144 * t147;
	t149 = pkin(7) + qJ(3);
	t153 = t144 * t149;
	t152 = t147 * t148;
	t140 = qJ(2) + pkin(11);
	t151 = cos(qJ(2));
	t150 = sin(qJ(2));
	t146 = cos(pkin(11));
	t142 = sin(pkin(11));
	t139 = cos(t140);
	t138 = sin(t140);
	t137 = pkin(6) - t140;
	t136 = pkin(6) + t140;
	t135 = -t142 * pkin(3) + qJ(4) * t146;
	t134 = pkin(3) * t146 + qJ(4) * t142 + pkin(2);
	t133 = cos(t136) + cos(t137);
	t132 = t138 * t152 + t143 * t139;
	t131 = t138 * t156 - t147 * t139;
	t1 = [-t131 * t145 + t143 * t157, t131 * t141 + t143 * t155, t147 * t138 + t143 * t133 / 0.2e1, (t147 * t134 + t135 * t156) * t151 + (-t134 * t156 + t147 * t135) * t150 + t143 * t153 + t147 * pkin(1) + 0; t132 * t145 - t141 * t154, -t132 * t141 - t145 * t154, t143 * t138 - t147 * t133 / 0.2e1, (t143 * t134 - t135 * t152) * t151 + (t134 * t152 + t143 * t135) * t150 - t147 * t153 + t143 * pkin(1) + 0; t138 * t155 + t148 * t141, -t138 * t157 + t148 * t145, -sin(t137) / 0.2e1 - sin(t136) / 0.2e1, t148 * t149 + qJ(1) + 0 + (t134 * t150 - t135 * t151) * t144; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:57:17
	% EndTime: 2020-11-04 20:57:17
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (108->39), mult. (98->55), div. (0->0), fcn. (118->18), ass. (0->30)
	t172 = qJ(2) + pkin(11);
	t168 = sin(t172);
	t175 = sin(pkin(6));
	t186 = t168 * t175;
	t174 = sin(pkin(10));
	t185 = t174 * t175;
	t178 = cos(pkin(6));
	t184 = t174 * t178;
	t177 = cos(pkin(10));
	t183 = t175 * t177;
	t182 = t177 * t178;
	t181 = cos(qJ(2));
	t180 = sin(qJ(2));
	t179 = pkin(8) + qJ(4);
	t176 = cos(pkin(11));
	t173 = sin(pkin(11));
	t171 = pkin(12) + qJ(5);
	t170 = cos(t172);
	t169 = cos(t171);
	t167 = sin(t171);
	t166 = pkin(6) - t172;
	t165 = pkin(6) + t172;
	t164 = cos(pkin(12)) * pkin(4) + pkin(3);
	t163 = sin(pkin(12)) * pkin(4) + qJ(3) + pkin(7);
	t162 = cos(t165) + cos(t166);
	t161 = -t173 * t164 + t179 * t176;
	t160 = t164 * t176 + t179 * t173 + pkin(2);
	t159 = t168 * t182 + t174 * t170;
	t158 = t168 * t184 - t177 * t170;
	t1 = [-t158 * t169 + t167 * t185, t158 * t167 + t169 * t185, t177 * t168 + t174 * t162 / 0.2e1, (t177 * t160 + t161 * t184) * t181 + (-t160 * t184 + t161 * t177) * t180 + t163 * t185 + t177 * pkin(1) + 0; t159 * t169 - t167 * t183, -t159 * t167 - t169 * t183, t174 * t168 - t177 * t162 / 0.2e1, (t160 * t174 - t161 * t182) * t181 + (t160 * t182 + t161 * t174) * t180 - t163 * t183 + t174 * pkin(1) + 0; t178 * t167 + t169 * t186, -t167 * t186 + t178 * t169, -sin(t166) / 0.2e1 - sin(t165) / 0.2e1, t163 * t178 + qJ(1) + 0 + (t160 * t180 - t161 * t181) * t175; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:57:17
	% EndTime: 2020-11-04 20:57:17
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (187->52), mult. (169->73), div. (0->0), fcn. (211->20), ass. (0->41)
	t210 = qJ(2) + pkin(11);
	t206 = sin(t210);
	t213 = sin(pkin(6));
	t226 = t206 * t213;
	t212 = sin(pkin(10));
	t225 = t212 * t213;
	t216 = cos(pkin(6));
	t224 = t212 * t216;
	t215 = cos(pkin(10));
	t223 = t213 * t215;
	t222 = t215 * t216;
	t221 = cos(qJ(2));
	t220 = cos(qJ(6));
	t219 = sin(qJ(2));
	t218 = sin(qJ(6));
	t217 = pkin(8) + qJ(4);
	t214 = cos(pkin(11));
	t211 = sin(pkin(11));
	t209 = pkin(12) + qJ(5);
	t208 = cos(t210);
	t207 = cos(t209);
	t205 = sin(t209);
	t204 = pkin(6) - t210;
	t203 = pkin(6) + t210;
	t202 = cos(pkin(12)) * pkin(4) + pkin(3);
	t201 = sin(pkin(12)) * pkin(4) + qJ(3) + pkin(7);
	t200 = cos(t203) + cos(t204);
	t199 = -sin(t204) / 0.2e1 - sin(t203) / 0.2e1;
	t198 = -t211 * t202 + t217 * t214;
	t197 = t202 * t214 + t217 * t211 + pkin(2);
	t196 = t206 * t222 + t212 * t208;
	t195 = t206 * t224 - t215 * t208;
	t194 = t216 * t205 + t207 * t226;
	t193 = t205 * t226 - t216 * t207;
	t192 = t215 * t206 + t212 * t200 / 0.2e1;
	t191 = t212 * t206 - t215 * t200 / 0.2e1;
	t190 = -t195 * t207 + t205 * t225;
	t189 = t196 * t207 - t205 * t223;
	t188 = t196 * t205 + t207 * t223;
	t187 = t195 * t205 + t207 * t225;
	t1 = [t190 * t220 + t192 * t218, -t190 * t218 + t192 * t220, -t187, t190 * pkin(5) - t187 * pkin(9) + (t215 * t197 + t198 * t224) * t221 + (-t197 * t224 + t198 * t215) * t219 + t201 * t225 + t215 * pkin(1) + 0; t189 * t220 + t191 * t218, -t189 * t218 + t191 * t220, t188, t189 * pkin(5) + t188 * pkin(9) + (t197 * t212 - t198 * t222) * t221 + (t197 * t222 + t198 * t212) * t219 - t201 * t223 + t212 * pkin(1) + 0; t194 * t220 + t199 * t218, -t194 * t218 + t199 * t220, t193, t194 * pkin(5) + t193 * pkin(9) + t201 * t216 + qJ(1) + 0 + (t197 * t219 - t198 * t221) * t213; 0, 0, 0, 1;];
	Tc_mdh = t1;
end