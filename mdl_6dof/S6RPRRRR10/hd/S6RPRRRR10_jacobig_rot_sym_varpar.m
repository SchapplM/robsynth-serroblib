% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRRR10
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:11
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RPRRRR10_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRR10_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_jacobig_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:11:29
	% EndTime: 2019-10-10 09:11:29
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:11:29
	% EndTime: 2019-10-10 09:11:29
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:11:29
	% EndTime: 2019-10-10 09:11:29
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:11:29
	% EndTime: 2019-10-10 09:11:30
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (5->5), mult. (15->13), div. (0->0), fcn. (23->8), ass. (0->11)
	t78 = sin(pkin(6));
	t80 = cos(pkin(7));
	t85 = t78 * t80;
	t79 = cos(pkin(13));
	t81 = cos(pkin(6));
	t84 = t79 * t81;
	t83 = cos(qJ(1));
	t82 = sin(qJ(1));
	t77 = sin(pkin(7));
	t76 = sin(pkin(13));
	t1 = [0, 0, -(-t83 * t76 - t82 * t84) * t77 + t82 * t85, 0, 0, 0; 0, 0, -(-t82 * t76 + t83 * t84) * t77 - t83 * t85, 0, 0, 0; 1, 0, -t78 * t79 * t77 + t81 * t80, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:11:30
	% EndTime: 2019-10-10 09:11:30
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (15->13), mult. (46->31), div. (0->0), fcn. (67->10), ass. (0->19)
	t121 = sin(pkin(6));
	t126 = sin(qJ(1));
	t134 = t121 * t126;
	t128 = cos(qJ(1));
	t133 = t121 * t128;
	t119 = sin(pkin(13));
	t132 = t126 * t119;
	t122 = cos(pkin(13));
	t131 = t126 * t122;
	t130 = t128 * t119;
	t129 = t128 * t122;
	t127 = cos(qJ(3));
	t125 = sin(qJ(3));
	t124 = cos(pkin(6));
	t123 = cos(pkin(7));
	t120 = sin(pkin(7));
	t118 = -t124 * t131 - t130;
	t117 = t124 * t129 - t132;
	t1 = [0, 0, -t118 * t120 + t123 * t134, (-t124 * t132 + t129) * t125 + (-t118 * t123 - t120 * t134) * t127, 0, 0; 0, 0, -t117 * t120 - t123 * t133, (t124 * t130 + t131) * t125 + (-t117 * t123 + t120 * t133) * t127, 0, 0; 1, 0, -t121 * t122 * t120 + t124 * t123, -t124 * t120 * t127 + (-t122 * t123 * t127 + t119 * t125) * t121, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:11:30
	% EndTime: 2019-10-10 09:11:30
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (25->13), mult. (77->31), div. (0->0), fcn. (111->10), ass. (0->22)
	t140 = sin(pkin(6));
	t145 = sin(qJ(1));
	t153 = t140 * t145;
	t147 = cos(qJ(1));
	t152 = t140 * t147;
	t138 = sin(pkin(13));
	t151 = t145 * t138;
	t141 = cos(pkin(13));
	t150 = t145 * t141;
	t149 = t147 * t138;
	t148 = t147 * t141;
	t146 = cos(qJ(3));
	t144 = sin(qJ(3));
	t143 = cos(pkin(6));
	t142 = cos(pkin(7));
	t139 = sin(pkin(7));
	t137 = -t143 * t150 - t149;
	t136 = t143 * t148 - t151;
	t135 = -t143 * t139 * t146 + (-t141 * t142 * t146 + t138 * t144) * t140;
	t134 = (-t143 * t151 + t148) * t144 + (-t137 * t142 - t139 * t153) * t146;
	t133 = (t143 * t149 + t150) * t144 + (-t136 * t142 + t139 * t152) * t146;
	t1 = [0, 0, -t137 * t139 + t142 * t153, t134, t134, 0; 0, 0, -t136 * t139 - t142 * t152, t133, t133, 0; 1, 0, -t140 * t141 * t139 + t143 * t142, t135, t135, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:11:30
	% EndTime: 2019-10-10 09:11:31
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (49->21), mult. (129->45), div. (0->0), fcn. (184->12), ass. (0->34)
	t192 = sin(pkin(7));
	t196 = cos(pkin(6));
	t210 = t192 * t196;
	t193 = sin(pkin(6));
	t198 = sin(qJ(1));
	t209 = t193 * t198;
	t200 = cos(qJ(1));
	t208 = t193 * t200;
	t194 = cos(pkin(13));
	t195 = cos(pkin(7));
	t207 = t194 * t195;
	t191 = sin(pkin(13));
	t206 = t198 * t191;
	t205 = t198 * t194;
	t204 = t200 * t191;
	t203 = t200 * t194;
	t184 = t196 * t203 - t206;
	t202 = -t184 * t195 + t192 * t208;
	t186 = -t196 * t205 - t204;
	t201 = t186 * t195 + t192 * t209;
	t199 = cos(qJ(3));
	t197 = sin(qJ(3));
	t190 = qJ(4) + qJ(5);
	t189 = cos(t190);
	t188 = sin(t190);
	t187 = -t196 * t206 + t203;
	t185 = t196 * t204 + t205;
	t183 = -t193 * t194 * t192 + t196 * t195;
	t182 = -t186 * t192 + t195 * t209;
	t181 = -t184 * t192 - t195 * t208;
	t180 = -t199 * t210 + (t191 * t197 - t199 * t207) * t193;
	t179 = t187 * t197 - t201 * t199;
	t178 = t185 * t197 + t202 * t199;
	t1 = [0, 0, t182, t179, t179, (t187 * t199 + t201 * t197) * t188 - t182 * t189; 0, 0, t181, t178, t178, (t185 * t199 - t202 * t197) * t188 - t181 * t189; 1, 0, t183, t180, t180, (t197 * t210 + (t191 * t199 + t197 * t207) * t193) * t188 - t183 * t189;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end