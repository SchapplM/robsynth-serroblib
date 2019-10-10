% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRRPRR11
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:10
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRPRR11_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR11_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_jacobig_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:10:18
	% EndTime: 2019-10-10 12:10:18
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:10:18
	% EndTime: 2019-10-10 12:10:18
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:10:18
	% EndTime: 2019-10-10 12:10:18
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t43 = sin(pkin(6));
	t1 = [0, sin(qJ(1)) * t43, 0, 0, 0, 0; 0, -cos(qJ(1)) * t43, 0, 0, 0, 0; 1, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:10:18
	% EndTime: 2019-10-10 12:10:18
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (4->4), mult. (9->8), div. (0->0), fcn. (17->6), ass. (0->8)
	t70 = cos(pkin(6));
	t73 = cos(qJ(2));
	t75 = t70 * t73;
	t74 = cos(qJ(1));
	t72 = sin(qJ(1));
	t71 = sin(qJ(2));
	t69 = sin(pkin(6));
	t1 = [0, t72 * t69, t74 * t71 + t72 * t75, 0, 0, 0; 0, -t74 * t69, t72 * t71 - t74 * t75, 0, 0, 0; 1, t70, -t69 * t73, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:10:19
	% EndTime: 2019-10-10 12:10:19
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (4->4), mult. (9->8), div. (0->0), fcn. (17->6), ass. (0->8)
	t90 = cos(pkin(6));
	t93 = cos(qJ(2));
	t95 = t90 * t93;
	t94 = cos(qJ(1));
	t92 = sin(qJ(1));
	t91 = sin(qJ(2));
	t89 = sin(pkin(6));
	t1 = [0, t92 * t89, t94 * t91 + t92 * t95, 0, 0, 0; 0, -t94 * t89, t92 * t91 - t94 * t95, 0, 0, 0; 1, t90, -t89 * t93, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:10:19
	% EndTime: 2019-10-10 12:10:19
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (6->6), mult. (16->8), div. (0->0), fcn. (29->6), ass. (0->11)
	t94 = sin(pkin(6));
	t98 = cos(qJ(2));
	t101 = t94 * t98;
	t95 = cos(pkin(6));
	t100 = t95 * t98;
	t99 = cos(qJ(1));
	t97 = sin(qJ(1));
	t96 = sin(qJ(2));
	t93 = t97 * t100 + t99 * t96;
	t92 = -t99 * t100 + t97 * t96;
	t1 = [0, t97 * t94, t93, 0, -t93, 0; 0, -t99 * t94, t92, 0, -t92, 0; 1, t95, -t101, 0, t101, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:10:19
	% EndTime: 2019-10-10 12:10:19
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (19->17), mult. (52->30), div. (0->0), fcn. (81->10), ass. (0->23)
	t177 = sin(pkin(6));
	t181 = sin(qJ(2));
	t194 = t177 * t181;
	t185 = cos(qJ(2));
	t193 = t177 * t185;
	t182 = sin(qJ(1));
	t192 = t182 * t177;
	t191 = t182 * t181;
	t190 = t182 * t185;
	t186 = cos(qJ(1));
	t189 = t186 * t177;
	t188 = t186 * t181;
	t187 = t186 * t185;
	t184 = cos(qJ(3));
	t183 = cos(qJ(5));
	t180 = sin(qJ(3));
	t179 = sin(qJ(5));
	t178 = cos(pkin(6));
	t176 = -t178 * t191 + t187;
	t175 = t178 * t190 + t188;
	t174 = t178 * t188 + t190;
	t173 = -t178 * t187 + t191;
	t1 = [0, t192, t175, 0, -t175, (t176 * t184 + t180 * t192) * t179 - (t176 * t180 - t184 * t192) * t183; 0, -t189, t173, 0, -t173, (t174 * t184 - t180 * t189) * t179 - (t174 * t180 + t184 * t189) * t183; 1, t178, -t193, 0, t193, (t178 * t180 + t184 * t194) * t179 - (-t178 * t184 + t180 * t194) * t183;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end