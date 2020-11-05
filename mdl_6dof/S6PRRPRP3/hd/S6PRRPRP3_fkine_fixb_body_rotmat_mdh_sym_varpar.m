% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRPRP3 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:09
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PRRPRP3_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP3_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRP3_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP3_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:09:08
	% EndTime: 2020-11-04 21:09:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:09:08
	% EndTime: 2020-11-04 21:09:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t101 = cos(pkin(10));
	t100 = sin(pkin(10));
	t1 = [t101, -t100, 0, 0; t100, t101, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:09:08
	% EndTime: 2020-11-04 21:09:08
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t102 = sin(pkin(10));
	t103 = sin(pkin(6));
	t111 = t102 * t103;
	t104 = cos(pkin(10));
	t110 = t104 * t103;
	t105 = cos(pkin(6));
	t106 = sin(qJ(2));
	t109 = t105 * t106;
	t107 = cos(qJ(2));
	t108 = t105 * t107;
	t1 = [-t102 * t109 + t104 * t107, -t102 * t108 - t104 * t106, t111, t104 * pkin(1) + pkin(7) * t111 + 0; t102 * t107 + t104 * t109, -t102 * t106 + t104 * t108, -t110, t102 * pkin(1) - pkin(7) * t110 + 0; t103 * t106, t103 * t107, t105, t105 * pkin(7) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:09:08
	% EndTime: 2020-11-04 21:09:08
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (29->27), mult. (64->48), div. (0->0), fcn. (85->8), ass. (0->18)
	t115 = sin(pkin(6));
	t128 = t115 * pkin(7);
	t114 = sin(pkin(10));
	t117 = cos(pkin(6));
	t127 = t114 * t117;
	t118 = sin(qJ(3));
	t126 = t115 * t118;
	t120 = cos(qJ(3));
	t125 = t115 * t120;
	t116 = cos(pkin(10));
	t124 = t116 * t117;
	t119 = sin(qJ(2));
	t123 = t117 * t119;
	t121 = cos(qJ(2));
	t122 = t117 * t121;
	t113 = t114 * t121 + t116 * t123;
	t112 = t114 * t123 - t116 * t121;
	t1 = [-t112 * t120 + t114 * t126, t112 * t118 + t114 * t125, t114 * t122 + t116 * t119, (t116 * pkin(2) + pkin(8) * t127) * t121 + (-pkin(2) * t127 + t116 * pkin(8)) * t119 + t114 * t128 + t116 * pkin(1) + 0; t113 * t120 - t116 * t126, -t113 * t118 - t116 * t125, t114 * t119 - t116 * t122, (t114 * pkin(2) - pkin(8) * t124) * t121 + (pkin(2) * t124 + t114 * pkin(8)) * t119 - t116 * t128 + t114 * pkin(1) + 0; t117 * t118 + t119 * t125, t117 * t120 - t119 * t126, -t115 * t121, t117 * pkin(7) + qJ(1) + 0 + (pkin(2) * t119 - pkin(8) * t121) * t115; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:09:08
	% EndTime: 2020-11-04 21:09:08
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (57->39), mult. (134->66), div. (0->0), fcn. (178->10), ass. (0->29)
	t141 = sin(pkin(6));
	t156 = t141 * pkin(7);
	t140 = sin(pkin(10));
	t144 = cos(pkin(6));
	t155 = t140 * t144;
	t145 = sin(qJ(3));
	t154 = t141 * t145;
	t147 = cos(qJ(3));
	t153 = t141 * t147;
	t148 = cos(qJ(2));
	t152 = t141 * t148;
	t143 = cos(pkin(10));
	t151 = t143 * t144;
	t146 = sin(qJ(2));
	t150 = t144 * t146;
	t149 = t144 * t148;
	t142 = cos(pkin(11));
	t139 = sin(pkin(11));
	t138 = t144 * t145 + t146 * t153;
	t137 = -t144 * t147 + t146 * t154;
	t136 = t140 * t149 + t143 * t146;
	t135 = t140 * t148 + t143 * t150;
	t134 = t140 * t146 - t143 * t149;
	t133 = t140 * t150 - t143 * t148;
	t132 = -t133 * t147 + t140 * t154;
	t131 = t135 * t147 - t143 * t154;
	t130 = t135 * t145 + t143 * t153;
	t129 = t133 * t145 + t140 * t153;
	t1 = [t132 * t142 + t136 * t139, -t132 * t139 + t136 * t142, -t129, t132 * pkin(3) - t129 * qJ(4) + (t143 * pkin(2) + pkin(8) * t155) * t148 + (-pkin(2) * t155 + t143 * pkin(8)) * t146 + t140 * t156 + t143 * pkin(1) + 0; t131 * t142 + t134 * t139, -t131 * t139 + t134 * t142, t130, t131 * pkin(3) + t130 * qJ(4) + (t140 * pkin(2) - pkin(8) * t151) * t148 + (pkin(2) * t151 + t140 * pkin(8)) * t146 - t143 * t156 + t140 * pkin(1) + 0; t138 * t142 - t139 * t152, -t138 * t139 - t142 * t152, t137, t138 * pkin(3) + t144 * pkin(7) + t137 * qJ(4) + qJ(1) + 0 + (pkin(2) * t146 - pkin(8) * t148) * t141; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:09:08
	% EndTime: 2020-11-04 21:09:08
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (90->36), mult. (149->58), div. (0->0), fcn. (183->12), ass. (0->33)
	t169 = sin(pkin(6));
	t171 = cos(pkin(6));
	t163 = sin(pkin(11)) * pkin(4) + pkin(8);
	t174 = sin(qJ(2));
	t176 = cos(qJ(2));
	t164 = cos(pkin(11)) * pkin(4) + pkin(3);
	t172 = qJ(4) + pkin(9);
	t173 = sin(qJ(3));
	t175 = cos(qJ(3));
	t180 = t164 * t175 + t172 * t173 + pkin(2);
	t178 = -t163 * t176 + t180 * t174;
	t179 = t164 * t173 - t172 * t175 + pkin(7);
	t191 = t179 * t169 - t178 * t171;
	t187 = t169 * t173;
	t186 = t169 * t175;
	t185 = t169 * t176;
	t184 = t171 * t174;
	t183 = t171 * t176;
	t182 = t174 * t175;
	t181 = t175 * t176;
	t177 = t163 * t174 + t180 * t176 + pkin(1);
	t170 = cos(pkin(10));
	t168 = sin(pkin(10));
	t167 = pkin(11) + qJ(5);
	t166 = cos(t167);
	t165 = sin(t167);
	t162 = t169 * t182 + t171 * t173;
	t161 = -t171 * t182 + t187;
	t160 = t168 * t183 + t170 * t174;
	t159 = t168 * t174 - t170 * t183;
	t158 = t168 * t161 + t170 * t181;
	t157 = -t170 * t161 + t168 * t181;
	t1 = [t158 * t166 + t160 * t165, -t158 * t165 + t160 * t166, -(t168 * t184 - t170 * t176) * t173 - t168 * t186, t191 * t168 + t177 * t170 + 0; t157 * t166 + t159 * t165, -t157 * t165 + t159 * t166, (t168 * t176 + t170 * t184) * t173 + t170 * t186, t177 * t168 - t191 * t170 + 0; t162 * t166 - t165 * t185, -t162 * t165 - t166 * t185, -t171 * t175 + t174 * t187, t178 * t169 + t179 * t171 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:09:08
	% EndTime: 2020-11-04 21:09:08
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (128->44), mult. (211->64), div. (0->0), fcn. (265->12), ass. (0->39)
	t210 = sin(pkin(6));
	t212 = cos(pkin(6));
	t204 = sin(pkin(11)) * pkin(4) + pkin(8);
	t215 = sin(qJ(2));
	t217 = cos(qJ(2));
	t205 = cos(pkin(11)) * pkin(4) + pkin(3);
	t213 = qJ(4) + pkin(9);
	t214 = sin(qJ(3));
	t216 = cos(qJ(3));
	t221 = t205 * t216 + t213 * t214 + pkin(2);
	t219 = -t204 * t217 + t221 * t215;
	t220 = t205 * t214 - t213 * t216 + pkin(7);
	t232 = t220 * t210 - t219 * t212;
	t228 = t210 * t214;
	t227 = t210 * t216;
	t226 = t210 * t217;
	t225 = t212 * t215;
	t224 = t212 * t217;
	t223 = t215 * t216;
	t222 = t216 * t217;
	t218 = t204 * t215 + t221 * t217 + pkin(1);
	t211 = cos(pkin(10));
	t209 = sin(pkin(10));
	t208 = pkin(11) + qJ(5);
	t207 = cos(t208);
	t206 = sin(t208);
	t203 = t210 * t223 + t212 * t214;
	t202 = -t212 * t223 + t228;
	t201 = t209 * t224 + t211 * t215;
	t200 = t209 * t215 - t211 * t224;
	t199 = t209 * t202 + t211 * t222;
	t198 = -t211 * t202 + t209 * t222;
	t197 = t203 * t207 - t206 * t226;
	t196 = t203 * t206 + t207 * t226;
	t195 = -t199 * t206 + t201 * t207;
	t194 = t199 * t207 + t201 * t206;
	t193 = -t198 * t206 + t200 * t207;
	t192 = t198 * t207 + t200 * t206;
	t1 = [t194, -(t209 * t225 - t211 * t217) * t214 - t209 * t227, -t195, t194 * pkin(5) - t195 * qJ(6) + t232 * t209 + t218 * t211 + 0; t192, (t209 * t217 + t211 * t225) * t214 + t211 * t227, -t193, t192 * pkin(5) - t193 * qJ(6) + t218 * t209 - t232 * t211 + 0; t197, -t212 * t216 + t215 * t228, t196, t197 * pkin(5) + t196 * qJ(6) + t219 * t210 + t220 * t212 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end