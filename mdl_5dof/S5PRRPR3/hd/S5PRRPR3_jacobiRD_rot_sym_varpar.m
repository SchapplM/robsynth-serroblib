% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PRRPR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5PRRPR3_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR3_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR3_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRPR3_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR3_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:20:25
	% EndTime: 2019-12-05 16:20:25
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:20:25
	% EndTime: 2019-12-05 16:20:25
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:20:25
	% EndTime: 2019-12-05 16:20:25
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (7->4), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->4)
	t35 = pkin(8) + qJ(2);
	t37 = qJD(2) * sin(t35);
	t36 = qJD(2) * cos(t35);
	t1 = [0, -t36, 0, 0, 0; 0, -t37, 0, 0, 0; 0, 0, 0, 0, 0; 0, t37, 0, 0, 0; 0, -t36, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:20:25
	% EndTime: 2019-12-05 16:20:26
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (29->10), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t41 = sin(qJ(3));
	t46 = qJD(2) * t41;
	t42 = cos(qJ(3));
	t45 = qJD(2) * t42;
	t44 = qJD(3) * t41;
	t43 = qJD(3) * t42;
	t40 = pkin(8) + qJ(2);
	t39 = cos(t40);
	t38 = sin(t40);
	t37 = t38 * t44 - t39 * t45;
	t36 = t38 * t43 + t39 * t46;
	t35 = t38 * t45 + t39 * t44;
	t34 = t38 * t46 - t39 * t43;
	t1 = [0, t37, t34, 0, 0; 0, -t35, -t36, 0, 0; 0, 0, -t44, 0, 0; 0, t36, t35, 0, 0; 0, t34, t37, 0, 0; 0, 0, -t43, 0, 0; 0, -qJD(2) * t38, 0, 0, 0; 0, qJD(2) * t39, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:20:25
	% EndTime: 2019-12-05 16:20:25
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (47->11), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->16)
	t56 = pkin(8) + qJ(2);
	t52 = sin(t56);
	t62 = qJD(2) * t52;
	t54 = cos(t56);
	t61 = qJD(2) * t54;
	t57 = qJ(3) + pkin(9);
	t55 = cos(t57);
	t60 = qJD(2) * t55;
	t53 = sin(t57);
	t59 = qJD(3) * t53;
	t58 = qJD(3) * t55;
	t51 = t52 * t59 - t54 * t60;
	t50 = t52 * t58 + t53 * t61;
	t49 = t52 * t60 + t54 * t59;
	t48 = t53 * t62 - t54 * t58;
	t1 = [0, t51, t48, 0, 0; 0, -t49, -t50, 0, 0; 0, 0, -t59, 0, 0; 0, t50, t49, 0, 0; 0, t48, t51, 0, 0; 0, 0, -t58, 0, 0; 0, -t62, 0, 0, 0; 0, t61, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:20:25
	% EndTime: 2019-12-05 16:20:26
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (115->15), mult. (54->12), div. (0->0), fcn. (54->4), ass. (0->16)
	t84 = qJ(3) + pkin(9) + qJ(5);
	t80 = sin(t84);
	t86 = qJD(3) + qJD(5);
	t90 = t86 * t80;
	t81 = cos(t84);
	t89 = t86 * t81;
	t85 = pkin(8) + qJ(2);
	t82 = sin(t85);
	t88 = qJD(2) * t82;
	t83 = cos(t85);
	t87 = qJD(2) * t83;
	t79 = -t81 * t87 + t82 * t90;
	t78 = t80 * t87 + t82 * t89;
	t77 = t81 * t88 + t83 * t90;
	t76 = t80 * t88 - t83 * t89;
	t1 = [0, t79, t76, 0, t76; 0, -t77, -t78, 0, -t78; 0, 0, -t90, 0, -t90; 0, t78, t77, 0, t77; 0, t76, t79, 0, t79; 0, 0, -t89, 0, -t89; 0, -t88, 0, 0, 0; 0, t87, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end