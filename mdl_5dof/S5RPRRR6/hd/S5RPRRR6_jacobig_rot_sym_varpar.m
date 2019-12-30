% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S5RPRRR6
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% Jg_rot [3x5]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 17:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S5RPRRR6_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR6_jacobig_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRR6_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR6_jacobig_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:43:40
	% EndTime: 2019-12-29 17:43:40
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:43:40
	% EndTime: 2019-12-29 17:43:40
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:43:34
	% EndTime: 2019-12-29 17:43:34
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:43:35
	% EndTime: 2019-12-29 17:43:35
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (3->2), mult. (0->0), div. (0->0), fcn. (2->2), ass. (0->2)
	t10 = qJ(1) + pkin(9);
	t1 = [0, 0, sin(t10), 0, 0; 0, 0, -cos(t10), 0, 0; 1, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:43:40
	% EndTime: 2019-12-29 17:43:40
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (6->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->4)
	t21 = qJ(1) + pkin(9);
	t20 = cos(t21);
	t19 = sin(t21);
	t1 = [0, 0, t19, t19, 0; 0, 0, -t20, -t20, 0; 1, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:43:35
	% EndTime: 2019-12-29 17:43:35
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (12->5), mult. (2->2), div. (0->0), fcn. (9->4), ass. (0->6)
	t86 = qJ(3) + qJ(4);
	t85 = qJ(1) + pkin(9);
	t84 = sin(t86);
	t83 = cos(t85);
	t82 = sin(t85);
	t1 = [0, 0, t82, t82, t83 * t84; 0, 0, -t83, -t83, t82 * t84; 1, 0, 0, 0, -cos(t86);];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,5);
end