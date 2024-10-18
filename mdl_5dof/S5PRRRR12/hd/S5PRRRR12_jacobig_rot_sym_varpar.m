% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S5PRRRR12
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
%   Siehe auch: S5PRRRR12_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d2,d3,d4,d5,theta1]';
% 
% Output:
% Jg_rot [3x5]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-28 18:09
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S5PRRRR12_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR12_jacobig_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRR12_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR12_jacobig_rot_sym_varpar: pkin has to be [11x1] (double)');
Jg_rot=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-28 18:09:12
	% EndTime: 2024-09-28 18:09:12
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-28 18:09:12
	% EndTime: 2024-09-28 18:09:12
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-28 18:09:12
	% EndTime: 2024-09-28 18:09:12
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t18 = sin(pkin(5));
	t1 = [0, sin(pkin(11)) * t18, 0, 0, 0; 0, -cos(pkin(11)) * t18, 0, 0, 0; 0, cos(pkin(5)), 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-28 18:09:12
	% EndTime: 2024-09-28 18:09:12
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (4->2), div. (0->0), fcn. (10->4), ass. (0->5)
	t40 = sin(pkin(5));
	t43 = cos(pkin(11)) * t40;
	t42 = cos(pkin(5));
	t39 = sin(pkin(11)) * t40;
	t1 = [0, t39, t39, 0, 0; 0, -t43, -t43, 0, 0; 0, t42, t42, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-28 18:09:12
	% EndTime: 2024-09-28 18:09:12
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (6->2), div. (0->0), fcn. (15->4), ass. (0->5)
	t52 = sin(pkin(5));
	t55 = cos(pkin(11)) * t52;
	t54 = cos(pkin(5));
	t51 = sin(pkin(11)) * t52;
	t1 = [0, t51, t51, t51, 0; 0, -t55, -t55, -t55, 0; 0, t54, t54, t54, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-28 18:09:12
	% EndTime: 2024-09-28 18:09:12
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (18->9), mult. (23->15), div. (0->0), fcn. (40->8), ass. (0->13)
	t150 = qJ(2) + qJ(3) + qJ(4);
	t149 = cos(t150);
	t156 = cos(pkin(5));
	t158 = t149 * t156;
	t153 = sin(pkin(5));
	t154 = cos(pkin(11));
	t157 = t154 * t153;
	t155 = cos(pkin(6));
	t152 = sin(pkin(6));
	t151 = sin(pkin(11));
	t148 = sin(t150);
	t147 = t151 * t153;
	t1 = [0, t147, t147, t147, t152 * t154 * t148 + (t152 * t158 + t153 * t155) * t151; 0, -t157, -t157, -t157, -t155 * t157 + (t148 * t151 - t154 * t158) * t152; 0, t156, t156, t156, -t153 * t149 * t152 + t156 * t155;];
	Jg_rot = t1;
end