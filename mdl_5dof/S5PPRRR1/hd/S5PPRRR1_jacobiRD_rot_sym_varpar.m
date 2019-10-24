% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PPRRR1
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
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:20
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5PPRRR1_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR1_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR1_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRRR1_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR1_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:20:29
	% EndTime: 2019-10-24 10:20:29
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:20:29
	% EndTime: 2019-10-24 10:20:29
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:20:29
	% EndTime: 2019-10-24 10:20:29
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:20:29
	% EndTime: 2019-10-24 10:20:29
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (10->5), mult. (10->8), div. (0->0), fcn. (10->4), ass. (0->6)
	t24 = qJD(3) * sin(pkin(8));
	t23 = qJD(3) * cos(pkin(8));
	t20 = pkin(9) + qJ(3);
	t19 = cos(t20);
	t18 = sin(t20);
	t1 = [0, 0, -t19 * t23, 0, 0; 0, 0, -t19 * t24, 0, 0; 0, 0, -qJD(3) * t18, 0, 0; 0, 0, t18 * t23, 0, 0; 0, 0, t18 * t24, 0, 0; 0, 0, -qJD(3) * t19, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:20:29
	% EndTime: 2019-10-24 10:20:29
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (44->10), mult. (20->8), div. (0->0), fcn. (20->4), ass. (0->13)
	t42 = pkin(9) + qJ(3) + qJ(4);
	t40 = sin(t42);
	t43 = qJD(3) + qJD(4);
	t51 = t43 * t40;
	t41 = cos(t42);
	t50 = t43 * t41;
	t49 = t43 * sin(pkin(8));
	t48 = t43 * cos(pkin(8));
	t47 = t41 * t49;
	t46 = t41 * t48;
	t39 = t40 * t48;
	t38 = t40 * t49;
	t1 = [0, 0, -t46, -t46, 0; 0, 0, -t47, -t47, 0; 0, 0, -t51, -t51, 0; 0, 0, t39, t39, 0; 0, 0, t38, t38, 0; 0, 0, -t50, -t50, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:20:30
	% EndTime: 2019-10-24 10:20:30
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (134->19), mult. (114->36), div. (0->0), fcn. (114->6), ass. (0->29)
	t232 = pkin(9) + qJ(3) + qJ(4);
	t230 = sin(t232);
	t233 = qJD(3) + qJD(4);
	t250 = t230 * t233;
	t236 = sin(qJ(5));
	t249 = t233 * t236;
	t237 = cos(qJ(5));
	t248 = t233 * t237;
	t234 = sin(pkin(8));
	t247 = t234 * t236;
	t246 = t234 * t237;
	t235 = cos(pkin(8));
	t245 = t235 * t236;
	t244 = t235 * t237;
	t243 = qJD(5) * t236;
	t242 = qJD(5) * t237;
	t241 = t234 * t250;
	t240 = t235 * t250;
	t231 = cos(t232);
	t239 = t230 * t242 + t231 * t249;
	t238 = t230 * t243 - t231 * t248;
	t229 = t233 * t231;
	t228 = -t230 * t248 - t231 * t243;
	t227 = t230 * t249 - t231 * t242;
	t226 = t238 * t235;
	t225 = t239 * t235;
	t224 = t238 * t234;
	t223 = t239 * t234;
	t1 = [0, 0, t226, t226, t236 * t240 + (-t231 * t244 - t247) * qJD(5); 0, 0, t224, t224, t236 * t241 + (-t231 * t246 + t245) * qJD(5); 0, 0, t228, t228, -t239; 0, 0, t225, t225, t237 * t240 + (t231 * t245 - t246) * qJD(5); 0, 0, t223, t223, t237 * t241 + (t231 * t247 + t244) * qJD(5); 0, 0, t227, t227, t238; 0, 0, -t240, -t240, 0; 0, 0, -t241, -t241, 0; 0, 0, t229, t229, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end