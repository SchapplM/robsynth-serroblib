% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRPRR7 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:12
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PRRPRR7_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRR7_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:12:36
	% EndTime: 2020-11-04 21:12:36
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:12:36
	% EndTime: 2020-11-04 21:12:36
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t112 = cos(pkin(11));
	t111 = sin(pkin(11));
	t1 = [t112, -t111, 0, 0; t111, t112, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:12:36
	% EndTime: 2020-11-04 21:12:36
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t113 = sin(pkin(11));
	t114 = sin(pkin(6));
	t122 = t113 * t114;
	t115 = cos(pkin(11));
	t121 = t115 * t114;
	t116 = cos(pkin(6));
	t117 = sin(qJ(2));
	t120 = t116 * t117;
	t118 = cos(qJ(2));
	t119 = t116 * t118;
	t1 = [-t113 * t120 + t115 * t118, -t113 * t119 - t115 * t117, t122, t115 * pkin(1) + pkin(7) * t122 + 0; t113 * t118 + t115 * t120, -t113 * t117 + t115 * t119, -t121, t113 * pkin(1) - pkin(7) * t121 + 0; t114 * t117, t114 * t118, t116, t116 * pkin(7) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:12:36
	% EndTime: 2020-11-04 21:12:36
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (29->27), mult. (64->48), div. (0->0), fcn. (85->8), ass. (0->18)
	t126 = sin(pkin(6));
	t139 = t126 * pkin(7);
	t125 = sin(pkin(11));
	t128 = cos(pkin(6));
	t138 = t125 * t128;
	t129 = sin(qJ(3));
	t137 = t126 * t129;
	t131 = cos(qJ(3));
	t136 = t126 * t131;
	t127 = cos(pkin(11));
	t135 = t127 * t128;
	t130 = sin(qJ(2));
	t134 = t128 * t130;
	t132 = cos(qJ(2));
	t133 = t128 * t132;
	t124 = t125 * t132 + t127 * t134;
	t123 = t125 * t134 - t127 * t132;
	t1 = [-t123 * t131 + t125 * t137, t123 * t129 + t125 * t136, t125 * t133 + t127 * t130, (t127 * pkin(2) + pkin(8) * t138) * t132 + (-pkin(2) * t138 + t127 * pkin(8)) * t130 + t125 * t139 + t127 * pkin(1) + 0; t124 * t131 - t127 * t137, -t124 * t129 - t127 * t136, t125 * t130 - t127 * t133, (t125 * pkin(2) - pkin(8) * t135) * t132 + (pkin(2) * t135 + t125 * pkin(8)) * t130 - t127 * t139 + t125 * pkin(1) + 0; t128 * t129 + t130 * t136, t128 * t131 - t130 * t137, -t126 * t132, t128 * pkin(7) + qJ(1) + 0 + (pkin(2) * t130 - pkin(8) * t132) * t126; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:12:36
	% EndTime: 2020-11-04 21:12:36
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (45->26), mult. (102->41), div. (0->0), fcn. (123->8), ass. (0->20)
	t143 = sin(pkin(6));
	t145 = cos(pkin(6));
	t147 = sin(qJ(2));
	t149 = cos(qJ(2));
	t146 = sin(qJ(3));
	t148 = cos(qJ(3));
	t153 = pkin(3) * t148 + qJ(4) * t146 + pkin(2);
	t151 = -pkin(8) * t149 + t153 * t147;
	t152 = t146 * pkin(3) - qJ(4) * t148 + pkin(7);
	t161 = t152 * t143 - t151 * t145;
	t157 = t145 * t147;
	t156 = t145 * t149;
	t155 = t146 * t143;
	t154 = t148 * t143;
	t150 = pkin(8) * t147 + t153 * t149 + pkin(1);
	t144 = cos(pkin(11));
	t142 = sin(pkin(11));
	t141 = t142 * t149 + t144 * t157;
	t140 = t142 * t157 - t144 * t149;
	t1 = [t142 * t156 + t144 * t147, t140 * t148 - t142 * t155, -t140 * t146 - t142 * t154, t161 * t142 + t150 * t144 + 0; t142 * t147 - t144 * t156, -t141 * t148 + t144 * t155, t141 * t146 + t144 * t154, t150 * t142 - t161 * t144 + 0; -t143 * t149, -t145 * t146 - t147 * t154, -t145 * t148 + t147 * t155, t151 * t143 + t152 * t145 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:12:36
	% EndTime: 2020-11-04 21:12:36
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (70->34), mult. (136->56), div. (0->0), fcn. (170->10), ass. (0->31)
	t169 = sin(pkin(6));
	t171 = cos(pkin(6));
	t174 = sin(qJ(2));
	t177 = cos(qJ(2));
	t178 = pkin(4) + pkin(8);
	t173 = sin(qJ(3));
	t176 = cos(qJ(3));
	t179 = pkin(3) + pkin(9);
	t183 = qJ(4) * t173 + t179 * t176 + pkin(2);
	t181 = t183 * t174 - t178 * t177;
	t182 = qJ(4) * t176 - t179 * t173 - pkin(7);
	t194 = t182 * t169 + t181 * t171;
	t191 = t169 * t173;
	t190 = t169 * t177;
	t189 = t171 * t174;
	t188 = t171 * t177;
	t187 = t173 * t174;
	t186 = t173 * t177;
	t185 = t176 * t169;
	t180 = t178 * t174 + t183 * t177 + pkin(1);
	t175 = cos(qJ(5));
	t172 = sin(qJ(5));
	t170 = cos(pkin(11));
	t168 = sin(pkin(11));
	t167 = t169 * t187 - t171 * t176;
	t166 = t171 * t187 + t185;
	t165 = t168 * t188 + t170 * t174;
	t164 = t168 * t174 - t170 * t188;
	t163 = -t168 * t166 + t170 * t186;
	t162 = t170 * t166 + t168 * t186;
	t1 = [t163 * t172 + t165 * t175, t163 * t175 - t165 * t172, -(t168 * t189 - t170 * t177) * t176 + t168 * t191, -t194 * t168 + t180 * t170 + 0; t162 * t172 + t164 * t175, t162 * t175 - t164 * t172, (t168 * t177 + t170 * t189) * t176 - t170 * t191, t180 * t168 + t194 * t170 + 0; t167 * t172 - t175 * t190, t167 * t175 + t172 * t190, t171 * t173 + t174 * t185, t181 * t169 - t182 * t171 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:12:36
	% EndTime: 2020-11-04 21:12:36
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (103->37), mult. (155->58), div. (0->0), fcn. (189->12), ass. (0->33)
	t207 = sin(pkin(6));
	t209 = cos(pkin(6));
	t201 = sin(qJ(5)) * pkin(5) + qJ(4);
	t204 = pkin(3) + pkin(9) + pkin(10);
	t210 = sin(qJ(3));
	t213 = cos(qJ(3));
	t218 = t201 * t213 - t204 * t210 - pkin(7);
	t211 = sin(qJ(2));
	t214 = cos(qJ(2));
	t219 = t201 * t210 + t204 * t213 + pkin(2);
	t220 = cos(qJ(5)) * pkin(5) + pkin(4) + pkin(8);
	t231 = -t219 * t211 + t220 * t214;
	t232 = -t218 * t207 + t231 * t209;
	t227 = t207 * t210;
	t226 = t207 * t214;
	t225 = t209 * t211;
	t224 = t209 * t214;
	t223 = t210 * t211;
	t222 = t210 * t214;
	t221 = t213 * t207;
	t216 = t220 * t211 + t219 * t214 + pkin(1);
	t208 = cos(pkin(11));
	t206 = sin(pkin(11));
	t205 = qJ(5) + qJ(6);
	t203 = cos(t205);
	t202 = sin(t205);
	t200 = t207 * t223 - t209 * t213;
	t199 = t209 * t223 + t221;
	t198 = t206 * t224 + t208 * t211;
	t197 = t206 * t211 - t208 * t224;
	t196 = -t206 * t199 + t208 * t222;
	t195 = t208 * t199 + t206 * t222;
	t1 = [t196 * t202 + t198 * t203, t196 * t203 - t198 * t202, -(t206 * t225 - t208 * t214) * t213 + t206 * t227, t232 * t206 + t216 * t208 + 0; t195 * t202 + t197 * t203, t195 * t203 - t197 * t202, (t206 * t214 + t208 * t225) * t213 - t208 * t227, t216 * t206 - t232 * t208 + 0; t200 * t202 - t203 * t226, t200 * t203 + t202 * t226, t209 * t210 + t211 * t221, -t231 * t207 - t218 * t209 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end