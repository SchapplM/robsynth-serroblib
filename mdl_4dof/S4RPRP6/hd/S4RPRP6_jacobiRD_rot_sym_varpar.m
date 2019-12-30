% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S4RPRP6
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
% 
% Output:
% JRD_rot [9x4]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 12:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S4RPRP6_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP6_jacobiRD_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP6_jacobiRD_rot_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RPRP6_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP6_jacobiRD_rot_sym_varpar: pkin has to be [5x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:58:02
	% EndTime: 2019-12-29 12:58:02
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:58:02
	% EndTime: 2019-12-29 12:58:02
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0; -t31, 0, 0, 0; 0, 0, 0, 0; t31, 0, 0, 0; -t30, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:58:02
	% EndTime: 2019-12-29 12:58:02
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (1->1), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t13 = qJD(1) * sin(qJ(1));
	t11 = qJD(1) * cos(qJ(1));
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; t11, 0, 0, 0; t13, 0, 0, 0; 0, 0, 0, 0; -t13, 0, 0, 0; t11, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:57:57
	% EndTime: 2019-12-29 12:57:57
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (11->9), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t35 = sin(qJ(1));
	t42 = qJD(1) * t35;
	t37 = cos(qJ(1));
	t41 = qJD(1) * t37;
	t34 = sin(qJ(3));
	t40 = qJD(3) * t34;
	t36 = cos(qJ(3));
	t39 = qJD(3) * t36;
	t38 = qJD(3) * t37;
	t33 = -t35 * t40 + t36 * t41;
	t32 = t34 * t41 + t35 * t39;
	t31 = t34 * t38 + t36 * t42;
	t30 = -t34 * t42 + t36 * t38;
	t1 = [t30, 0, t33, 0; t32, 0, t31, 0; 0, 0, -t39, 0; -t31, 0, -t32, 0; t33, 0, t30, 0; 0, 0, t40, 0; -t41, 0, 0, 0; -t42, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:58:02
	% EndTime: 2019-12-29 12:58:02
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (11->9), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t45 = sin(qJ(1));
	t52 = qJD(1) * t45;
	t47 = cos(qJ(1));
	t51 = qJD(1) * t47;
	t44 = sin(qJ(3));
	t50 = qJD(3) * t44;
	t46 = cos(qJ(3));
	t49 = qJD(3) * t46;
	t48 = qJD(3) * t47;
	t43 = -t45 * t50 + t46 * t51;
	t42 = t44 * t51 + t45 * t49;
	t41 = t44 * t48 + t46 * t52;
	t40 = -t44 * t52 + t46 * t48;
	t1 = [t40, 0, t43, 0; t42, 0, t41, 0; 0, 0, -t49, 0; -t41, 0, -t42, 0; t43, 0, t40, 0; 0, 0, t50, 0; -t51, 0, 0, 0; -t52, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,4);
end