% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPRPR3
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
%   Siehe auch: S5RPRPR3_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:21
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RPRPR3_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR3_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
JRD_rot=NaN(9,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:21:20
	% EndTime: 2022-01-23 09:21:20
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:21:20
	% EndTime: 2022-01-23 09:21:20
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
	% StartTime: 2022-01-23 09:21:20
	% EndTime: 2022-01-23 09:21:20
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
	% StartTime: 2022-01-23 09:21:20
	% EndTime: 2022-01-23 09:21:20
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
	% StartTime: 2022-01-23 09:21:20
	% EndTime: 2022-01-23 09:21:20
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
	% StartTime: 2022-01-23 09:21:20
	% EndTime: 2022-01-23 09:21:21
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (173->20), mult. (132->28), div. (0->0), fcn. (132->6), ass. (0->24)
	t199 = qJ(1) + pkin(8) + qJ(3);
	t197 = sin(t199);
	t204 = cos(qJ(5));
	t215 = t197 * t204;
	t198 = cos(t199);
	t203 = sin(qJ(5));
	t214 = t198 * t203;
	t200 = qJD(1) + qJD(3);
	t201 = sin(pkin(9));
	t213 = t200 * t201;
	t202 = cos(pkin(9));
	t212 = t202 * t203;
	t211 = t202 * t204;
	t210 = qJD(5) * t203;
	t209 = qJD(5) * t204;
	t208 = t197 * t213;
	t207 = t198 * t213;
	t206 = -t197 * t203 - t198 * t211;
	t205 = t197 * t212 + t198 * t204;
	t194 = t205 * qJD(5) + t206 * t200;
	t193 = (t198 * t212 - t215) * t200 + (t197 * t211 - t214) * qJD(5);
	t192 = -t200 * t214 - t197 * t209 + (t198 * t210 + t200 * t215) * t202;
	t191 = t206 * qJD(5) + t205 * t200;
	t1 = [t194, 0, t194, 0, t191; -t192, 0, -t192, 0, -t193; 0, 0, 0, 0, -t201 * t209; t193, 0, t193, 0, t192; t191, 0, t191, 0, t194; 0, 0, 0, 0, t201 * t210; -t207, 0, -t207, 0, 0; -t208, 0, -t208, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
end