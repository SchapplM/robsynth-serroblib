% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 15:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RPPRP3_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP3_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRP3_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_jacobiRD_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:59:26
	% EndTime: 2019-12-29 15:59:26
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:59:31
	% EndTime: 2019-12-29 15:59:31
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0; -t31, 0, 0, 0, 0; 0, 0, 0, 0, 0; t31, 0, 0, 0, 0; -t30, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:59:26
	% EndTime: 2019-12-29 15:59:26
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (7->4), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->4)
	t35 = qJ(1) + pkin(7);
	t37 = qJD(1) * sin(t35);
	t36 = qJD(1) * cos(t35);
	t1 = [-t36, 0, 0, 0, 0; -t37, 0, 0, 0, 0; 0, 0, 0, 0, 0; t37, 0, 0, 0, 0; -t36, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:59:32
	% EndTime: 2019-12-29 15:59:32
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (5->2), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->4)
	t16 = qJ(1) + pkin(7);
	t17 = qJD(1) * sin(t16);
	t14 = qJD(1) * cos(t16);
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; t14, 0, 0, 0, 0; t17, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t17, 0, 0, 0, 0; t14, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:59:27
	% EndTime: 2019-12-29 15:59:27
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (29->10), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t42 = sin(qJ(4));
	t47 = qJD(1) * t42;
	t43 = cos(qJ(4));
	t46 = qJD(1) * t43;
	t45 = qJD(4) * t42;
	t44 = qJD(4) * t43;
	t41 = qJ(1) + pkin(7);
	t40 = cos(t41);
	t39 = sin(t41);
	t38 = -t39 * t45 + t40 * t46;
	t37 = t39 * t44 + t40 * t47;
	t36 = t39 * t46 + t40 * t45;
	t35 = -t39 * t47 + t40 * t44;
	t1 = [t35, 0, 0, t38, 0; t37, 0, 0, t36, 0; 0, 0, 0, -t44, 0; -t36, 0, 0, -t37, 0; t38, 0, 0, t35, 0; 0, 0, 0, t45, 0; -qJD(1) * t40, 0, 0, 0, 0; -qJD(1) * t39, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:59:32
	% EndTime: 2019-12-29 15:59:32
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (29->10), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t52 = sin(qJ(4));
	t57 = qJD(1) * t52;
	t53 = cos(qJ(4));
	t56 = qJD(1) * t53;
	t55 = qJD(4) * t52;
	t54 = qJD(4) * t53;
	t51 = qJ(1) + pkin(7);
	t50 = cos(t51);
	t49 = sin(t51);
	t48 = -t49 * t55 + t50 * t56;
	t47 = t49 * t54 + t50 * t57;
	t46 = t49 * t56 + t50 * t55;
	t45 = -t49 * t57 + t50 * t54;
	t1 = [t45, 0, 0, t48, 0; t47, 0, 0, t46, 0; 0, 0, 0, -t54, 0; -t46, 0, 0, -t47, 0; t48, 0, 0, t45, 0; 0, 0, 0, t55, 0; -qJD(1) * t50, 0, 0, 0, 0; -qJD(1) * t49, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end