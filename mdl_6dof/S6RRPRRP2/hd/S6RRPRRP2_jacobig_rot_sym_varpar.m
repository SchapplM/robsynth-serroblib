% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:31
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRPRRP2_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP2_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRP2_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP2_jacobig_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:31:51
	% EndTime: 2019-10-10 10:31:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:31:51
	% EndTime: 2019-10-10 10:31:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:31:51
	% EndTime: 2019-10-10 10:31:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, sin(qJ(1)), 0, 0, 0, 0; 0, -cos(qJ(1)), 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:31:51
	% EndTime: 2019-10-10 10:31:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, sin(qJ(1)), 0, 0, 0, 0; 0, -cos(qJ(1)), 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:31:51
	% EndTime: 2019-10-10 10:31:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t22 = cos(qJ(1));
	t21 = sin(qJ(1));
	t1 = [0, t21, 0, t21, 0, 0; 0, -t22, 0, -t22, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:31:51
	% EndTime: 2019-10-10 10:31:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->4), mult. (2->2), div. (0->0), fcn. (9->4), ass. (0->5)
	t91 = cos(qJ(1));
	t90 = sin(qJ(1));
	t89 = qJ(2) + pkin(10) + qJ(4);
	t88 = sin(t89);
	t1 = [0, t90, 0, t90, t91 * t88, 0; 0, -t91, 0, -t91, t90 * t88, 0; 1, 0, 0, 0, -cos(t89), 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:31:51
	% EndTime: 2019-10-10 10:31:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->4), mult. (2->2), div. (0->0), fcn. (9->4), ass. (0->5)
	t106 = cos(qJ(1));
	t105 = sin(qJ(1));
	t104 = qJ(2) + pkin(10) + qJ(4);
	t103 = sin(t104);
	t1 = [0, t105, 0, t105, t106 * t103, 0; 0, -t106, 0, -t106, t105 * t103, 0; 1, 0, 0, 0, -cos(t104), 0;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end