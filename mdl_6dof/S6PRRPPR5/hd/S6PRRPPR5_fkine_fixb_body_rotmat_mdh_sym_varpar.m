% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRPPR5 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:08
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PRRPPR5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPPR5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:08:00
	% EndTime: 2020-11-04 21:08:00
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:08:00
	% EndTime: 2020-11-04 21:08:00
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t114 = cos(pkin(10));
	t113 = sin(pkin(10));
	t1 = [t114, -t113, 0, 0; t113, t114, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:08:00
	% EndTime: 2020-11-04 21:08:00
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t115 = sin(pkin(10));
	t116 = sin(pkin(6));
	t124 = t115 * t116;
	t117 = cos(pkin(10));
	t123 = t117 * t116;
	t118 = cos(pkin(6));
	t119 = sin(qJ(2));
	t122 = t118 * t119;
	t120 = cos(qJ(2));
	t121 = t118 * t120;
	t1 = [-t115 * t122 + t117 * t120, -t115 * t121 - t117 * t119, t124, t117 * pkin(1) + pkin(7) * t124 + 0; t115 * t120 + t117 * t122, -t115 * t119 + t117 * t121, -t123, t115 * pkin(1) - pkin(7) * t123 + 0; t116 * t119, t116 * t120, t118, t118 * pkin(7) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:08:00
	% EndTime: 2020-11-04 21:08:01
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (29->27), mult. (64->48), div. (0->0), fcn. (85->8), ass. (0->18)
	t128 = sin(pkin(6));
	t141 = t128 * pkin(7);
	t127 = sin(pkin(10));
	t130 = cos(pkin(6));
	t140 = t127 * t130;
	t131 = sin(qJ(3));
	t139 = t128 * t131;
	t133 = cos(qJ(3));
	t138 = t128 * t133;
	t129 = cos(pkin(10));
	t137 = t129 * t130;
	t132 = sin(qJ(2));
	t136 = t130 * t132;
	t134 = cos(qJ(2));
	t135 = t130 * t134;
	t126 = -t127 * t136 + t129 * t134;
	t125 = t127 * t134 + t129 * t136;
	t1 = [t126 * t133 + t127 * t139, -t126 * t131 + t127 * t138, t127 * t135 + t129 * t132, (t129 * pkin(2) + pkin(8) * t140) * t134 + (-pkin(2) * t140 + t129 * pkin(8)) * t132 + t127 * t141 + t129 * pkin(1) + 0; t125 * t133 - t129 * t139, -t125 * t131 - t129 * t138, t127 * t132 - t129 * t135, (t127 * pkin(2) - pkin(8) * t137) * t134 + (pkin(2) * t137 + t127 * pkin(8)) * t132 - t129 * t141 + t127 * pkin(1) + 0; t130 * t131 + t132 * t138, t130 * t133 - t132 * t139, -t128 * t134, t130 * pkin(7) + qJ(1) + 0 + (pkin(2) * t132 - pkin(8) * t134) * t128; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:08:01
	% EndTime: 2020-11-04 21:08:01
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (45->26), mult. (102->41), div. (0->0), fcn. (123->8), ass. (0->20)
	t145 = sin(pkin(6));
	t147 = cos(pkin(6));
	t149 = sin(qJ(2));
	t151 = cos(qJ(2));
	t148 = sin(qJ(3));
	t150 = cos(qJ(3));
	t155 = pkin(3) * t150 + qJ(4) * t148 + pkin(2);
	t153 = -pkin(8) * t151 + t155 * t149;
	t154 = t148 * pkin(3) - qJ(4) * t150 + pkin(7);
	t163 = t154 * t145 - t153 * t147;
	t159 = t147 * t149;
	t158 = t147 * t151;
	t157 = t148 * t145;
	t156 = t150 * t145;
	t152 = pkin(8) * t149 + t155 * t151 + pkin(1);
	t146 = cos(pkin(10));
	t144 = sin(pkin(10));
	t143 = t144 * t159 - t146 * t151;
	t142 = t144 * t151 + t146 * t159;
	t1 = [t144 * t158 + t146 * t149, t143 * t150 - t144 * t157, -t143 * t148 - t144 * t156, t163 * t144 + t152 * t146 + 0; t144 * t149 - t146 * t158, -t142 * t150 + t146 * t157, t142 * t148 + t146 * t156, t152 * t144 - t163 * t146 + 0; -t145 * t151, -t147 * t148 - t149 * t156, -t147 * t150 + t149 * t157, t153 * t145 + t154 * t147 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:08:01
	% EndTime: 2020-11-04 21:08:01
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (70->35), mult. (140->57), div. (0->0), fcn. (174->10), ass. (0->32)
	t171 = sin(pkin(6));
	t174 = cos(pkin(6));
	t177 = sin(qJ(2));
	t179 = cos(qJ(2));
	t180 = pkin(4) + pkin(8);
	t175 = qJ(5) + pkin(3);
	t176 = sin(qJ(3));
	t178 = cos(qJ(3));
	t184 = qJ(4) * t176 + t175 * t178 + pkin(2);
	t182 = t184 * t177 - t180 * t179;
	t183 = qJ(4) * t178 - t175 * t176 - pkin(7);
	t197 = t183 * t171 + t182 * t174;
	t169 = sin(pkin(11));
	t194 = t169 * t176;
	t170 = sin(pkin(10));
	t193 = t170 * t177;
	t192 = t171 * t176;
	t191 = t171 * t179;
	t172 = cos(pkin(11));
	t190 = t172 * t176;
	t173 = cos(pkin(10));
	t189 = t173 * t177;
	t188 = t174 * t177;
	t187 = t174 * t179;
	t186 = t178 * t171;
	t181 = t180 * t177 + t184 * t179 + pkin(1);
	t168 = -t174 * t178 + t177 * t192;
	t167 = t170 * t188 - t173 * t179;
	t166 = t170 * t179 + t173 * t188;
	t165 = t169 * t187 + t172 * t186;
	t164 = t169 * t186 - t172 * t187;
	t1 = [-t170 * t164 - t167 * t194 + t172 * t189, -t170 * t165 - t167 * t190 - t169 * t189, -t167 * t178 + t170 * t192, -t197 * t170 + t181 * t173 + 0; t173 * t164 + t166 * t194 + t172 * t193, t173 * t165 + t166 * t190 - t169 * t193, t166 * t178 - t173 * t192, t181 * t170 + t197 * t173 + 0; t168 * t169 - t172 * t191, t168 * t172 + t169 * t191, t174 * t176 + t177 * t186, t182 * t171 - t183 * t174 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:08:01
	% EndTime: 2020-11-04 21:08:01
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (103->37), mult. (149->58), div. (0->0), fcn. (183->12), ass. (0->33)
	t211 = sin(pkin(6));
	t213 = cos(pkin(6));
	t204 = cos(pkin(11)) * pkin(5) + pkin(4) + pkin(8);
	t215 = sin(qJ(2));
	t217 = cos(qJ(2));
	t205 = sin(pkin(11)) * pkin(5) + qJ(4);
	t208 = qJ(5) + pkin(3) + pkin(9);
	t214 = sin(qJ(3));
	t216 = cos(qJ(3));
	t221 = t205 * t214 + t208 * t216 + pkin(2);
	t219 = -t204 * t217 + t221 * t215;
	t220 = t205 * t216 - t208 * t214 - pkin(7);
	t232 = t220 * t211 + t219 * t213;
	t228 = t211 * t214;
	t227 = t211 * t216;
	t226 = t211 * t217;
	t225 = t213 * t215;
	t224 = t213 * t217;
	t223 = t214 * t215;
	t222 = t214 * t217;
	t218 = t204 * t215 + t221 * t217 + pkin(1);
	t212 = cos(pkin(10));
	t210 = sin(pkin(10));
	t209 = pkin(11) + qJ(6);
	t207 = cos(t209);
	t206 = sin(t209);
	t203 = t211 * t223 - t213 * t216;
	t202 = t213 * t223 + t227;
	t201 = t210 * t224 + t212 * t215;
	t200 = t210 * t215 - t212 * t224;
	t199 = -t210 * t202 + t212 * t222;
	t198 = t212 * t202 + t210 * t222;
	t1 = [t199 * t206 + t201 * t207, t199 * t207 - t201 * t206, -(t210 * t225 - t212 * t217) * t216 + t210 * t228, -t232 * t210 + t218 * t212 + 0; t198 * t206 + t200 * t207, t198 * t207 - t200 * t206, (t210 * t217 + t212 * t225) * t216 - t212 * t228, t218 * t210 + t232 * t212 + 0; t203 * t206 - t207 * t226, t203 * t207 + t206 * t226, t213 * t214 + t215 * t227, t219 * t211 - t220 * t213 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end