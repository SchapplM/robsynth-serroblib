% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRPP9 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:37
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRRRPP9_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP9_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPP9_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP9_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:37:32
	% EndTime: 2020-11-04 22:37:32
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:37:32
	% EndTime: 2020-11-04 22:37:32
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t106 = cos(qJ(1));
	t105 = sin(qJ(1));
	t1 = [t106, -t105, 0, 0; t105, t106, 0, 0; 0, 0, 1, pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:37:32
	% EndTime: 2020-11-04 22:37:32
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t107 = sin(pkin(6));
	t110 = sin(qJ(1));
	t118 = t110 * t107;
	t109 = sin(qJ(2));
	t117 = t110 * t109;
	t111 = cos(qJ(2));
	t116 = t110 * t111;
	t112 = cos(qJ(1));
	t115 = t112 * t107;
	t114 = t112 * t109;
	t113 = t112 * t111;
	t108 = cos(pkin(6));
	t1 = [-t108 * t117 + t113, -t108 * t116 - t114, t118, t112 * pkin(1) + pkin(8) * t118 + 0; t108 * t114 + t116, t108 * t113 - t117, -t115, t110 * pkin(1) - pkin(8) * t115 + 0; t107 * t109, t107 * t111, t108, t108 * pkin(8) + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:37:32
	% EndTime: 2020-11-04 22:37:32
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->23), mult. (60->37), div. (0->0), fcn. (81->8), ass. (0->20)
	t122 = sin(pkin(6));
	t127 = cos(qJ(3));
	t137 = t122 * t127;
	t124 = sin(qJ(3));
	t136 = t124 * t122;
	t125 = sin(qJ(2));
	t135 = t125 * t127;
	t126 = sin(qJ(1));
	t134 = t126 * t125;
	t128 = cos(qJ(2));
	t133 = t126 * t128;
	t129 = cos(qJ(1));
	t132 = t129 * t125;
	t131 = t129 * t128;
	t130 = pkin(2) * t125 - pkin(9) * t128;
	t123 = cos(pkin(6));
	t121 = t128 * pkin(2) + t125 * pkin(9) + pkin(1);
	t120 = -t123 * t132 - t133;
	t119 = t122 * pkin(8) - t130 * t123;
	t1 = [(-t123 * t135 + t136) * t126 + t127 * t131, (t123 * t134 - t131) * t124 + t126 * t137, t123 * t133 + t132, t119 * t126 + t121 * t129 + 0; -t120 * t127 - t129 * t136, t120 * t124 - t129 * t137, -t123 * t131 + t134, -t119 * t129 + t121 * t126 + 0; t122 * t135 + t123 * t124, t123 * t127 - t125 * t136, -t122 * t128, t123 * pkin(8) + t130 * t122 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:37:32
	% EndTime: 2020-11-04 22:37:32
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (57->34), mult. (112->60), div. (0->0), fcn. (146->10), ass. (0->29)
	t145 = sin(pkin(6));
	t152 = cos(qJ(3));
	t166 = t145 * t152;
	t147 = sin(qJ(4));
	t153 = cos(qJ(2));
	t165 = t147 * t153;
	t148 = sin(qJ(3));
	t164 = t148 * t145;
	t149 = sin(qJ(2));
	t163 = t149 * t152;
	t150 = sin(qJ(1));
	t162 = t150 * t153;
	t151 = cos(qJ(4));
	t161 = t151 * t153;
	t160 = t152 * t153;
	t154 = cos(qJ(1));
	t159 = t154 * t149;
	t158 = t154 * t153;
	t143 = t152 * pkin(3) + t148 * pkin(10) + pkin(2);
	t157 = t153 * pkin(9) - t143 * t149;
	t156 = t148 * pkin(3) - t152 * pkin(10) + pkin(8);
	t146 = cos(pkin(6));
	t141 = t146 * t163 - t164;
	t155 = t141 * t147 + t146 * t161;
	t142 = t147 * t160 - t151 * t149;
	t140 = t145 * t163 + t146 * t148;
	t139 = t149 * pkin(9) + t143 * t153 + pkin(1);
	t138 = t145 * t156 + t157 * t146;
	t1 = [(-t141 * t150 + t152 * t158) * t151 + (t146 * t162 + t159) * t147, -t154 * t142 + t155 * t150, (-t150 * t146 * t149 + t158) * t148 - t150 * t166, t138 * t150 + t139 * t154 + 0; (t141 * t151 - t146 * t165) * t154 + t150 * (t147 * t149 + t151 * t160), -t150 * t142 - t155 * t154, (t146 * t159 + t162) * t148 + t154 * t166, -t138 * t154 + t139 * t150 + 0; t140 * t151 - t145 * t165, -t140 * t147 - t145 * t161, -t146 * t152 + t149 * t164, -t157 * t145 + t156 * t146 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:37:32
	% EndTime: 2020-11-04 22:37:33
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (83->38), mult. (140->64), div. (0->0), fcn. (174->10), ass. (0->31)
	t179 = sin(qJ(4));
	t183 = cos(qJ(4));
	t173 = pkin(4) * t183 + qJ(5) * t179 + pkin(3);
	t180 = sin(qJ(3));
	t184 = cos(qJ(3));
	t200 = -t184 * pkin(10) + t173 * t180 + pkin(8);
	t177 = sin(pkin(6));
	t198 = t177 * t184;
	t185 = cos(qJ(2));
	t197 = t179 * t185;
	t196 = t180 * t177;
	t181 = sin(qJ(2));
	t195 = t181 * t184;
	t182 = sin(qJ(1));
	t194 = t182 * t185;
	t193 = t183 * t185;
	t192 = t184 * t185;
	t186 = cos(qJ(1));
	t191 = t186 * t181;
	t190 = t186 * t185;
	t169 = t180 * pkin(10) + t173 * t184 + pkin(2);
	t174 = -t179 * pkin(4) + qJ(5) * t183 - pkin(9);
	t188 = t169 * t181 + t174 * t185;
	t178 = cos(pkin(6));
	t171 = t178 * t195 - t196;
	t187 = t171 * t179 + t178 * t193;
	t172 = t179 * t192 - t183 * t181;
	t170 = t177 * t195 + t178 * t180;
	t168 = t169 * t185 - t174 * t181 + pkin(1);
	t167 = -t200 * t177 + t188 * t178;
	t1 = [(-t182 * t178 * t181 + t190) * t180 - t182 * t198, (t171 * t182 - t184 * t190) * t183 - (t178 * t194 + t191) * t179, t186 * t172 - t187 * t182, -t167 * t182 + t168 * t186 + 0; (t178 * t191 + t194) * t180 + t186 * t198, (-t171 * t183 + t178 * t197) * t186 - t182 * (t179 * t181 + t183 * t192), t182 * t172 + t187 * t186, t167 * t186 + t168 * t182 + 0; -t178 * t184 + t181 * t196, -t170 * t183 + t177 * t197, t170 * t179 + t177 * t193, t188 * t177 + t200 * t178 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:37:33
	% EndTime: 2020-11-04 22:37:33
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (104->40), mult. (142->64), div. (0->0), fcn. (176->10), ass. (0->33)
	t213 = qJ(6) + pkin(4);
	t214 = sin(qJ(4));
	t218 = cos(qJ(4));
	t207 = qJ(5) * t214 + t213 * t218 + pkin(3);
	t215 = sin(qJ(3));
	t219 = cos(qJ(3));
	t222 = pkin(5) + pkin(10);
	t203 = t207 * t219 + t222 * t215 + pkin(2);
	t216 = sin(qJ(2));
	t220 = cos(qJ(2));
	t224 = qJ(5) * t218 - t213 * t214 - pkin(9);
	t240 = t203 * t216 + t224 * t220;
	t238 = t207 * t215 - t222 * t219 + pkin(8);
	t211 = sin(pkin(6));
	t235 = t211 * t219;
	t234 = t214 * t220;
	t233 = t215 * t211;
	t232 = t216 * t219;
	t217 = sin(qJ(1));
	t231 = t217 * t220;
	t230 = t218 * t220;
	t229 = t219 * t220;
	t221 = cos(qJ(1));
	t228 = t221 * t216;
	t227 = t221 * t220;
	t212 = cos(pkin(6));
	t205 = t212 * t232 - t233;
	t223 = t205 * t214 + t212 * t230;
	t206 = t214 * t229 - t218 * t216;
	t204 = t211 * t232 + t212 * t215;
	t202 = t203 * t220 - t224 * t216 + pkin(1);
	t201 = -t211 * t238 + t240 * t212;
	t1 = [(-t217 * t212 * t216 + t227) * t215 - t217 * t235, t221 * t206 - t223 * t217, (-t205 * t217 + t219 * t227) * t218 + (t212 * t231 + t228) * t214, -t201 * t217 + t202 * t221 + 0; (t212 * t228 + t231) * t215 + t221 * t235, t217 * t206 + t223 * t221, (t205 * t218 - t212 * t234) * t221 + t217 * (t214 * t216 + t218 * t229), t201 * t221 + t202 * t217 + 0; -t212 * t219 + t216 * t233, t204 * t214 + t211 * t230, t204 * t218 - t211 * t234, t240 * t211 + t238 * t212 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end