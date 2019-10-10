% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RPPRRR4
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:06
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RPPRRR4_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR4_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRR4_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR4_jacobig_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:06:31
	% EndTime: 2019-10-10 00:06:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:06:31
	% EndTime: 2019-10-10 00:06:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:06:31
	% EndTime: 2019-10-10 00:06:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:06:31
	% EndTime: 2019-10-10 00:06:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:06:31
	% EndTime: 2019-10-10 00:06:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (8->4), ass. (0->5)
	t20 = cos(qJ(1));
	t19 = sin(qJ(1));
	t18 = cos(pkin(10));
	t17 = sin(pkin(10));
	t1 = [0, 0, 0, t20 * t17 - t19 * t18, 0, 0; 0, 0, 0, t19 * t17 + t20 * t18, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:06:32
	% EndTime: 2019-10-10 00:06:32
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (4->3), mult. (10->6), div. (0->0), fcn. (19->6), ass. (0->8)
	t75 = cos(qJ(1));
	t74 = sin(qJ(1));
	t73 = sin(qJ(4));
	t72 = cos(pkin(10));
	t71 = sin(pkin(10));
	t70 = t74 * t71 + t75 * t72;
	t69 = t75 * t71 - t74 * t72;
	t1 = [0, 0, 0, t69, t70 * t73, 0; 0, 0, 0, t70, -t69 * t73, 0; 1, 0, 0, 0, cos(qJ(4)), 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:06:32
	% EndTime: 2019-10-10 00:06:32
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (6->4), mult. (16->6), div. (0->0), fcn. (30->6), ass. (0->11)
	t99 = cos(qJ(1));
	t98 = cos(qJ(4));
	t97 = sin(qJ(1));
	t96 = sin(qJ(4));
	t95 = cos(pkin(10));
	t94 = sin(pkin(10));
	t93 = t97 * t94 + t99 * t95;
	t92 = t99 * t94 - t97 * t95;
	t91 = t93 * t96;
	t90 = t92 * t96;
	t1 = [0, 0, 0, t92, t91, t91; 0, 0, 0, t93, -t90, -t90; 1, 0, 0, 0, t98, t98;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end