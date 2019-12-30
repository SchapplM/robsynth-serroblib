% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S4RRPR3
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% Jg_rot [3x4]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 13:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S4RRPR3_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR3_jacobig_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RRPR3_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR3_jacobig_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:31:12
	% EndTime: 2019-12-29 13:31:12
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:31:12
	% EndTime: 2019-12-29 13:31:12
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 1, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:31:07
	% EndTime: 2019-12-29 13:31:07
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 1, 1, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:31:12
	% EndTime: 2019-12-29 13:31:12
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 1, 1, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:31:12
	% EndTime: 2019-12-29 13:31:12
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (5->2), mult. (0->0), div. (0->0), fcn. (2->2), ass. (0->2)
	t24 = qJ(1) + qJ(2) + pkin(7);
	t1 = [0, 0, 0, sin(t24); 0, 0, 0, -cos(t24); 1, 1, 0, 0;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,4);
end