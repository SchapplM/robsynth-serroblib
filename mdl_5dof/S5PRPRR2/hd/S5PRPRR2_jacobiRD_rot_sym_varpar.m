% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PRPRR2
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:25
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5PRPRR2_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR2_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR2_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPRR2_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR2_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:25:57
	% EndTime: 2019-10-24 10:25:57
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:25:57
	% EndTime: 2019-10-24 10:25:57
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:25:57
	% EndTime: 2019-10-24 10:25:57
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (4->4), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->5)
	t20 = qJD(2) * sin(qJ(2));
	t19 = qJD(2) * cos(qJ(2));
	t16 = cos(pkin(8));
	t15 = sin(pkin(8));
	t1 = [0, -t16 * t19, 0, 0, 0; 0, -t15 * t19, 0, 0, 0; 0, -t20, 0, 0, 0; 0, t16 * t20, 0, 0, 0; 0, t15 * t20, 0, 0, 0; 0, -t19, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:25:57
	% EndTime: 2019-10-24 10:25:57
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (10->5), mult. (10->8), div. (0->0), fcn. (10->4), ass. (0->6)
	t26 = qJD(2) * sin(pkin(8));
	t25 = qJD(2) * cos(pkin(8));
	t22 = qJ(2) + pkin(9);
	t21 = cos(t22);
	t20 = sin(t22);
	t1 = [0, -t21 * t25, 0, 0, 0; 0, -t21 * t26, 0, 0, 0; 0, -qJD(2) * t20, 0, 0, 0; 0, t20 * t25, 0, 0, 0; 0, t20 * t26, 0, 0, 0; 0, -qJD(2) * t21, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:25:57
	% EndTime: 2019-10-24 10:25:57
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (44->10), mult. (20->8), div. (0->0), fcn. (20->4), ass. (0->13)
	t41 = qJ(2) + pkin(9) + qJ(4);
	t39 = sin(t41);
	t42 = qJD(2) + qJD(4);
	t50 = t42 * t39;
	t40 = cos(t41);
	t49 = t42 * t40;
	t48 = t42 * sin(pkin(8));
	t47 = t42 * cos(pkin(8));
	t46 = t40 * t48;
	t45 = t40 * t47;
	t38 = t39 * t47;
	t37 = t39 * t48;
	t1 = [0, -t45, 0, -t45, 0; 0, -t46, 0, -t46, 0; 0, -t50, 0, -t50, 0; 0, t38, 0, t38, 0; 0, t37, 0, t37, 0; 0, -t49, 0, -t49, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:25:57
	% EndTime: 2019-10-24 10:25:57
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (134->19), mult. (114->36), div. (0->0), fcn. (114->6), ass. (0->29)
	t231 = qJ(2) + pkin(9) + qJ(4);
	t229 = sin(t231);
	t232 = qJD(2) + qJD(4);
	t249 = t229 * t232;
	t235 = sin(qJ(5));
	t248 = t232 * t235;
	t236 = cos(qJ(5));
	t247 = t232 * t236;
	t233 = sin(pkin(8));
	t246 = t233 * t235;
	t245 = t233 * t236;
	t234 = cos(pkin(8));
	t244 = t234 * t235;
	t243 = t234 * t236;
	t242 = qJD(5) * t235;
	t241 = qJD(5) * t236;
	t240 = t233 * t249;
	t239 = t234 * t249;
	t230 = cos(t231);
	t238 = t229 * t241 + t230 * t248;
	t237 = t229 * t242 - t230 * t247;
	t228 = t232 * t230;
	t227 = -t229 * t247 - t230 * t242;
	t226 = t229 * t248 - t230 * t241;
	t225 = t237 * t234;
	t224 = t238 * t234;
	t223 = t237 * t233;
	t222 = t238 * t233;
	t1 = [0, t225, 0, t225, t235 * t239 + (-t230 * t243 - t246) * qJD(5); 0, t223, 0, t223, t235 * t240 + (-t230 * t245 + t244) * qJD(5); 0, t227, 0, t227, -t238; 0, t224, 0, t224, t236 * t239 + (t230 * t244 - t245) * qJD(5); 0, t222, 0, t222, t236 * t240 + (t230 * t246 + t243) * qJD(5); 0, t226, 0, t226, t237; 0, -t239, 0, -t239, 0; 0, -t240, 0, -t240, 0; 0, t228, 0, t228, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end