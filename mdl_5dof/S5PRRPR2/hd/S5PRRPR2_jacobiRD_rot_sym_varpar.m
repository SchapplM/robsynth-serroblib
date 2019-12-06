% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PRRPR2
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
% Datum: 2019-12-05 16:18
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5PRRPR2_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR2_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR2_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRPR2_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR2_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:18:04
	% EndTime: 2019-12-05 16:18:04
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:18:04
	% EndTime: 2019-12-05 16:18:04
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:18:04
	% EndTime: 2019-12-05 16:18:04
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
	% StartTime: 2019-12-05 16:18:04
	% EndTime: 2019-12-05 16:18:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (30->8), mult. (8->2), div. (0->0), fcn. (8->2), ass. (0->5)
	t49 = pkin(8) + qJ(2) + qJ(3);
	t50 = qJD(2) + qJD(3);
	t51 = t50 * cos(t49);
	t46 = t50 * sin(t49);
	t1 = [0, -t51, -t51, 0, 0; 0, -t46, -t46, 0, 0; 0, 0, 0, 0, 0; 0, t46, t46, 0, 0; 0, -t51, -t51, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:18:04
	% EndTime: 2019-12-05 16:18:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (42->8), mult. (20->8), div. (0->0), fcn. (20->4), ass. (0->13)
	t56 = pkin(8) + qJ(2) + qJ(3);
	t54 = sin(t56);
	t57 = qJD(2) + qJD(3);
	t64 = t57 * t54;
	t63 = t57 * sin(pkin(9));
	t62 = t57 * cos(pkin(9));
	t61 = t54 * t62;
	t55 = cos(t56);
	t60 = t55 * t62;
	t53 = t57 * t55;
	t52 = t55 * t63;
	t51 = t54 * t63;
	t1 = [0, -t60, -t60, 0, 0; 0, -t61, -t61, 0, 0; 0, 0, 0, 0, 0; 0, t52, t52, 0, 0; 0, t51, t51, 0, 0; 0, 0, 0, 0, 0; 0, -t64, -t64, 0, 0; 0, t53, t53, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:18:04
	% EndTime: 2019-12-05 16:18:05
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (173->20), mult. (132->28), div. (0->0), fcn. (132->6), ass. (0->24)
	t200 = pkin(8) + qJ(2) + qJ(3);
	t198 = sin(t200);
	t205 = cos(qJ(5));
	t216 = t198 * t205;
	t199 = cos(t200);
	t204 = sin(qJ(5));
	t215 = t199 * t204;
	t201 = qJD(2) + qJD(3);
	t202 = sin(pkin(9));
	t214 = t201 * t202;
	t203 = cos(pkin(9));
	t213 = t203 * t204;
	t212 = t203 * t205;
	t211 = qJD(5) * t204;
	t210 = qJD(5) * t205;
	t209 = t198 * t214;
	t208 = t199 * t214;
	t207 = -t198 * t204 - t199 * t212;
	t206 = t198 * t213 + t199 * t205;
	t195 = t206 * qJD(5) + t207 * t201;
	t194 = (t199 * t213 - t216) * t201 + (t198 * t212 - t215) * qJD(5);
	t193 = -t201 * t215 - t198 * t210 + (t199 * t211 + t201 * t216) * t203;
	t192 = t207 * qJD(5) + t206 * t201;
	t1 = [0, t195, t195, 0, t192; 0, -t193, -t193, 0, -t194; 0, 0, 0, 0, -t202 * t210; 0, t194, t194, 0, t193; 0, t192, t192, 0, t195; 0, 0, 0, 0, t202 * t211; 0, -t208, -t208, 0, 0; 0, -t209, -t209, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end