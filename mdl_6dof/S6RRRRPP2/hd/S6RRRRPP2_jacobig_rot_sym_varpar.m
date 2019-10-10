% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRRRPP2
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:22
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRRPP2_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP2_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPP2_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP2_jacobig_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:22:15
	% EndTime: 2019-10-10 12:22:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:22:15
	% EndTime: 2019-10-10 12:22:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:22:15
	% EndTime: 2019-10-10 12:22:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, sin(qJ(1)), 0, 0, 0, 0; 0, -cos(qJ(1)), 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:22:15
	% EndTime: 2019-10-10 12:22:15
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t19 = cos(qJ(1));
	t18 = sin(qJ(1));
	t1 = [0, t18, t18, 0, 0, 0; 0, -t19, -t19, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:22:16
	% EndTime: 2019-10-10 12:22:16
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (6->4), mult. (2->2), div. (0->0), fcn. (9->4), ass. (0->5)
	t88 = cos(qJ(1));
	t87 = sin(qJ(1));
	t86 = qJ(2) + qJ(3);
	t85 = sin(t86);
	t1 = [0, t87, t87, t88 * t85, 0, 0; 0, -t88, -t88, t87 * t85, 0, 0; 1, 0, 0, -cos(t86), 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:22:16
	% EndTime: 2019-10-10 12:22:16
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (6->4), mult. (2->2), div. (0->0), fcn. (9->4), ass. (0->5)
	t103 = cos(qJ(1));
	t102 = sin(qJ(1));
	t101 = qJ(2) + qJ(3);
	t100 = sin(t101);
	t1 = [0, t102, t102, t103 * t100, 0, 0; 0, -t103, -t103, t102 * t100, 0, 0; 1, 0, 0, -cos(t101), 0, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:22:16
	% EndTime: 2019-10-10 12:22:16
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (6->4), mult. (2->2), div. (0->0), fcn. (9->4), ass. (0->5)
	t89 = cos(qJ(1));
	t88 = sin(qJ(1));
	t87 = qJ(2) + qJ(3);
	t86 = sin(t87);
	t1 = [0, t88, t88, t89 * t86, 0, 0; 0, -t89, -t89, t88 * t86, 0, 0; 1, 0, 0, -cos(t87), 0, 0;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end