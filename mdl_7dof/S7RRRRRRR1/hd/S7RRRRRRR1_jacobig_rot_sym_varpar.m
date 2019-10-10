% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
% 
% Output:
% Jg_rot [3x7]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 17:10
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S7RRRRRRR1_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),uint8(0),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_jacobig_rot_sym_varpar: qJ has to be [7x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S7RRRRRRR1_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_jacobig_rot_sym_varpar: pkin has to be [4x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 17:10:04
	% EndTime: 2019-10-10 17:10:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 17:10:04
	% EndTime: 2019-10-10 17:10:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 17:10:04
	% EndTime: 2019-10-10 17:10:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, sin(qJ(1)), 0, 0, 0, 0, 0; 0, -cos(qJ(1)), 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 17:10:04
	% EndTime: 2019-10-10 17:10:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (2->2), div. (0->0), fcn. (7->4), ass. (0->4)
	t52 = cos(qJ(1));
	t51 = sin(qJ(1));
	t50 = sin(qJ(2));
	t1 = [0, t51, -t52 * t50, 0, 0, 0, 0; 0, -t52, -t51 * t50, 0, 0, 0, 0; 1, 0, cos(qJ(2)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 17:10:04
	% EndTime: 2019-10-10 17:10:04
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (6->6), mult. (9->8), div. (0->0), fcn. (19->6), ass. (0->8)
	t76 = sin(qJ(3));
	t80 = cos(qJ(2));
	t82 = t76 * t80;
	t81 = cos(qJ(1));
	t79 = cos(qJ(3));
	t78 = sin(qJ(1));
	t77 = sin(qJ(2));
	t1 = [0, t78, -t81 * t77, -t78 * t79 - t81 * t82, 0, 0, 0; 0, -t81, -t78 * t77, -t78 * t82 + t81 * t79, 0, 0, 0; 1, 0, t80, -t77 * t76, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 17:10:05
	% EndTime: 2019-10-10 17:10:05
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (11->11), mult. (24->19), div. (0->0), fcn. (42->8), ass. (0->14)
	t128 = sin(qJ(2));
	t129 = sin(qJ(1));
	t138 = t129 * t128;
	t132 = cos(qJ(2));
	t137 = t129 * t132;
	t127 = sin(qJ(3));
	t133 = cos(qJ(1));
	t136 = t133 * t127;
	t135 = t133 * t128;
	t131 = cos(qJ(3));
	t134 = t133 * t131;
	t130 = cos(qJ(4));
	t126 = sin(qJ(4));
	t1 = [0, t129, -t135, -t129 * t131 - t132 * t136, (-t129 * t127 + t132 * t134) * t126 - t130 * t135, 0, 0; 0, -t133, -t138, -t127 * t137 + t134, (t131 * t137 + t136) * t126 - t130 * t138, 0, 0; 1, 0, t132, -t128 * t127, t128 * t131 * t126 + t132 * t130, 0, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 17:10:05
	% EndTime: 2019-10-10 17:10:05
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (21->17), mult. (52->31), div. (0->0), fcn. (83->10), ass. (0->22)
	t175 = sin(qJ(3));
	t176 = sin(qJ(2));
	t189 = t176 * t175;
	t180 = cos(qJ(3));
	t188 = t176 * t180;
	t177 = sin(qJ(1));
	t187 = t177 * t176;
	t181 = cos(qJ(2));
	t186 = t177 * t181;
	t182 = cos(qJ(1));
	t185 = t182 * t175;
	t184 = t182 * t176;
	t183 = t182 * t180;
	t179 = cos(qJ(4));
	t178 = cos(qJ(5));
	t174 = sin(qJ(4));
	t173 = sin(qJ(5));
	t172 = -t177 * t175 + t181 * t183;
	t171 = -t177 * t180 - t181 * t185;
	t170 = t180 * t186 + t185;
	t169 = -t175 * t186 + t183;
	t1 = [0, t177, -t184, t171, t172 * t174 - t179 * t184, (t172 * t179 + t174 * t184) * t173 - t171 * t178, 0; 0, -t182, -t187, t169, t170 * t174 - t179 * t187, (t170 * t179 + t174 * t187) * t173 - t169 * t178, 0; 1, 0, t181, -t189, t174 * t188 + t181 * t179, (-t181 * t174 + t179 * t188) * t173 + t178 * t189, 0;];
	Jg_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobig_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 17:10:07
	% EndTime: 2019-10-10 17:10:07
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (39->23), mult. (101->43), div. (0->0), fcn. (153->12), ass. (0->30)
	t254 = sin(qJ(3));
	t255 = sin(qJ(2));
	t269 = t255 * t254;
	t260 = cos(qJ(3));
	t268 = t255 * t260;
	t256 = sin(qJ(1));
	t267 = t256 * t255;
	t261 = cos(qJ(2));
	t266 = t256 * t261;
	t262 = cos(qJ(1));
	t265 = t262 * t254;
	t264 = t262 * t255;
	t263 = t262 * t260;
	t259 = cos(qJ(4));
	t258 = cos(qJ(5));
	t257 = cos(qJ(6));
	t253 = sin(qJ(4));
	t252 = sin(qJ(5));
	t251 = sin(qJ(6));
	t250 = -t256 * t254 + t261 * t263;
	t249 = -t256 * t260 - t261 * t265;
	t248 = t260 * t266 + t265;
	t247 = -t254 * t266 + t263;
	t246 = -t261 * t253 + t259 * t268;
	t245 = t253 * t268 + t261 * t259;
	t244 = t250 * t259 + t253 * t264;
	t243 = t250 * t253 - t259 * t264;
	t242 = t248 * t259 + t253 * t267;
	t241 = t248 * t253 - t259 * t267;
	t1 = [0, t256, -t264, t249, t243, t244 * t252 - t249 * t258, -(t244 * t258 + t249 * t252) * t251 + t243 * t257; 0, -t262, -t267, t247, t241, t242 * t252 - t247 * t258, -(t242 * t258 + t247 * t252) * t251 + t241 * t257; 1, 0, t261, -t269, t245, t246 * t252 + t258 * t269, -(t246 * t258 - t252 * t269) * t251 + t245 * t257;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,7);
end