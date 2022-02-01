% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPRPR2
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
%   Siehe auch: S5RPRPR2_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:19
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RPRPR2_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR2_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR2_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
JRD_rot=NaN(9,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:19:33
	% EndTime: 2022-01-23 09:19:33
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:19:33
	% EndTime: 2022-01-23 09:19:33
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0; -t31, 0, 0, 0, 0; 0, 0, 0, 0, 0; t31, 0, 0, 0, 0; -t30, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:19:33
	% EndTime: 2022-01-23 09:19:33
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (7->4), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->4)
	t35 = qJ(1) + pkin(8);
	t37 = qJD(1) * sin(t35);
	t36 = qJD(1) * cos(t35);
	t1 = [-t36, 0, 0, 0, 0; -t37, 0, 0, 0, 0; 0, 0, 0, 0, 0; t37, 0, 0, 0, 0; -t36, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:19:33
	% EndTime: 2022-01-23 09:19:33
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (30->8), mult. (8->2), div. (0->0), fcn. (8->2), ass. (0->5)
	t48 = qJ(1) + pkin(8) + qJ(3);
	t49 = qJD(1) + qJD(3);
	t50 = t49 * cos(t48);
	t45 = t49 * sin(t48);
	t1 = [-t50, 0, -t50, 0, 0; -t45, 0, -t45, 0, 0; 0, 0, 0, 0, 0; t45, 0, t45, 0, 0; -t50, 0, -t50, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:19:33
	% EndTime: 2022-01-23 09:19:33
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (42->8), mult. (20->8), div. (0->0), fcn. (20->4), ass. (0->13)
	t55 = qJ(1) + pkin(8) + qJ(3);
	t53 = sin(t55);
	t56 = qJD(1) + qJD(3);
	t63 = t56 * t53;
	t62 = t56 * sin(pkin(9));
	t61 = t56 * cos(pkin(9));
	t60 = t53 * t61;
	t54 = cos(t55);
	t59 = t54 * t61;
	t52 = t56 * t54;
	t51 = t54 * t62;
	t50 = t53 * t62;
	t1 = [-t59, 0, -t59, 0, 0; -t60, 0, -t60, 0, 0; 0, 0, 0, 0, 0; t51, 0, t51, 0, 0; t50, 0, t50, 0, 0; 0, 0, 0, 0, 0; -t63, 0, -t63, 0, 0; t52, 0, t52, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:19:33
	% EndTime: 2022-01-23 09:19:33
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (114->14), mult. (54->12), div. (0->0), fcn. (54->4), ass. (0->16)
	t77 = qJ(1) + pkin(8) + qJ(3);
	t73 = sin(t77);
	t79 = qJD(1) + qJD(3);
	t82 = t79 * t73;
	t74 = cos(t77);
	t72 = t79 * t74;
	t78 = pkin(9) + qJ(5);
	t75 = sin(t78);
	t81 = qJD(5) * t75;
	t76 = cos(t78);
	t80 = qJD(5) * t76;
	t71 = -t76 * t72 + t73 * t81;
	t70 = t75 * t72 + t73 * t80;
	t69 = t74 * t81 + t76 * t82;
	t68 = -t74 * t80 + t75 * t82;
	t1 = [t71, 0, t71, 0, t68; -t69, 0, -t69, 0, -t70; 0, 0, 0, 0, -t81; t70, 0, t70, 0, t69; t68, 0, t68, 0, t71; 0, 0, 0, 0, -t80; -t82, 0, -t82, 0, 0; t72, 0, t72, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
end