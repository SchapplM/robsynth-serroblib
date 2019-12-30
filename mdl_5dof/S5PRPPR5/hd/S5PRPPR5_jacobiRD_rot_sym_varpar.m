% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PRPPR5
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta4]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 15:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5PRPPR5_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR5_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR5_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPPR5_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR5_jacobiRD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:31:25
	% EndTime: 2019-12-29 15:31:25
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:31:30
	% EndTime: 2019-12-29 15:31:30
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:31:20
	% EndTime: 2019-12-29 15:31:20
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (4->4), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->5)
	t20 = qJD(2) * sin(qJ(2));
	t19 = qJD(2) * cos(qJ(2));
	t16 = cos(pkin(7));
	t15 = sin(pkin(7));
	t1 = [0, -t16 * t19, 0, 0, 0; 0, -t15 * t19, 0, 0, 0; 0, -t20, 0, 0, 0; 0, t16 * t20, 0, 0, 0; 0, t15 * t20, 0, 0, 0; 0, -t19, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:31:25
	% EndTime: 2019-12-29 15:31:25
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (5->5), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->5)
	t104 = qJD(2) * sin(qJ(2));
	t103 = qJD(2) * cos(qJ(2));
	t100 = cos(pkin(7));
	t99 = sin(pkin(7));
	t1 = [0, -t100 * t103, 0, 0, 0; 0, -t99 * t103, 0, 0, 0; 0, -t104, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, -t100 * t104, 0, 0, 0; 0, -t99 * t104, 0, 0, 0; 0, t103, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:31:25
	% EndTime: 2019-12-29 15:31:25
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (6->4), mult. (32->10), div. (0->0), fcn. (32->6), ass. (0->9)
	t28 = sin(pkin(8));
	t30 = cos(pkin(8));
	t32 = sin(qJ(2));
	t33 = cos(qJ(2));
	t35 = (t28 * t32 + t30 * t33) * qJD(2);
	t34 = (t28 * t33 - t30 * t32) * qJD(2);
	t31 = cos(pkin(7));
	t29 = sin(pkin(7));
	t1 = [0, -t31 * t35, 0, 0, 0; 0, -t29 * t35, 0, 0, 0; 0, t34, 0, 0, 0; 0, t31 * t34, 0, 0, 0; 0, t29 * t34, 0, 0, 0; 0, t35, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:31:26
	% EndTime: 2019-12-29 15:31:26
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (43->21), mult. (164->48), div. (0->0), fcn. (176->8), ass. (0->23)
	t239 = sin(pkin(8));
	t241 = cos(pkin(8));
	t244 = sin(qJ(2));
	t246 = cos(qJ(2));
	t235 = t244 * t239 + t246 * t241;
	t253 = qJD(2) * t235;
	t243 = sin(qJ(5));
	t250 = qJD(5) * t243;
	t245 = cos(qJ(5));
	t249 = qJD(5) * t245;
	t248 = t246 * t239 - t244 * t241;
	t234 = t248 * qJD(2);
	t242 = cos(pkin(7));
	t240 = sin(pkin(7));
	t232 = t235 * t242;
	t231 = t248 * t242;
	t230 = t235 * t240;
	t229 = t248 * t240;
	t228 = t242 * t234;
	t227 = t242 * t253;
	t226 = t240 * t234;
	t225 = t240 * t253;
	t1 = [0, -t227 * t245 - t231 * t250, 0, 0, -t228 * t243 + (-t232 * t245 + t240 * t243) * qJD(5); 0, -t225 * t245 - t229 * t250, 0, 0, -t226 * t243 + (-t230 * t245 - t242 * t243) * qJD(5); 0, t234 * t245 - t235 * t250, 0, 0, -t243 * t253 + t248 * t249; 0, t227 * t243 - t231 * t249, 0, 0, -t228 * t245 + (t232 * t243 + t240 * t245) * qJD(5); 0, t225 * t243 - t229 * t249, 0, 0, -t226 * t245 + (t230 * t243 - t242 * t245) * qJD(5); 0, -t234 * t243 - t235 * t249, 0, 0, -t245 * t253 - t248 * t250; 0, -t228, 0, 0, 0; 0, -t226, 0, 0, 0; 0, -t253, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end