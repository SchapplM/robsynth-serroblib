% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PRPRP6
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5PRPRP6_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP6_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP6_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPRP6_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRPRP6_jacobiRD_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:41:38
	% EndTime: 2019-12-05 15:41:38
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:41:38
	% EndTime: 2019-12-05 15:41:38
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:41:38
	% EndTime: 2019-12-05 15:41:38
	% DurationCPUTime: 0.04s
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
	% StartTime: 2019-12-05 15:41:38
	% EndTime: 2019-12-05 15:41:38
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (2->2), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->5)
	t103 = qJD(2) * sin(qJ(2));
	t102 = qJD(2) * cos(qJ(2));
	t99 = cos(pkin(7));
	t98 = sin(pkin(7));
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, t99 * t102, 0, 0, 0; 0, t98 * t102, 0, 0, 0; 0, t103, 0, 0, 0; 0, -t99 * t103, 0, 0, 0; 0, -t98 * t103, 0, 0, 0; 0, t102, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:41:38
	% EndTime: 2019-12-05 15:41:39
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (19->17), mult. (77->36), div. (0->0), fcn. (77->6), ass. (0->17)
	t148 = sin(qJ(4));
	t149 = sin(qJ(2));
	t161 = t148 * t149;
	t150 = cos(qJ(4));
	t160 = t149 * t150;
	t159 = qJD(2) * t149;
	t151 = cos(qJ(2));
	t158 = qJD(2) * t151;
	t157 = qJD(4) * t149;
	t156 = qJD(4) * t151;
	t146 = sin(pkin(7));
	t155 = t146 * t158;
	t147 = cos(pkin(7));
	t154 = t147 * t158;
	t153 = t148 * t156 + t150 * t159;
	t152 = -t148 * t159 + t150 * t156;
	t1 = [0, t152 * t147, 0, t150 * t154 + (-t146 * t150 - t147 * t161) * qJD(4), 0; 0, t152 * t146, 0, t150 * t155 + (-t146 * t161 + t147 * t150) * qJD(4), 0; 0, t148 * t158 + t150 * t157, 0, t153, 0; 0, -t153 * t147, 0, -t148 * t154 + (t146 * t148 - t147 * t160) * qJD(4), 0; 0, -t153 * t146, 0, -t148 * t155 + (-t146 * t160 - t147 * t148) * qJD(4), 0; 0, -t148 * t157 + t150 * t158, 0, t152, 0; 0, -t154, 0, 0, 0; 0, -t155, 0, 0, 0; 0, -t159, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:41:39
	% EndTime: 2019-12-05 15:41:39
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (19->17), mult. (77->36), div. (0->0), fcn. (77->6), ass. (0->17)
	t201 = sin(qJ(4));
	t202 = sin(qJ(2));
	t214 = t201 * t202;
	t203 = cos(qJ(4));
	t213 = t202 * t203;
	t212 = qJD(2) * t202;
	t204 = cos(qJ(2));
	t211 = qJD(2) * t204;
	t210 = qJD(4) * t202;
	t209 = qJD(4) * t204;
	t199 = sin(pkin(7));
	t208 = t199 * t211;
	t200 = cos(pkin(7));
	t207 = t200 * t211;
	t206 = t201 * t209 + t203 * t212;
	t205 = t201 * t212 - t203 * t209;
	t1 = [0, -t205 * t200, 0, t203 * t207 + (-t199 * t203 - t200 * t214) * qJD(4), 0; 0, -t205 * t199, 0, t203 * t208 + (-t199 * t214 + t200 * t203) * qJD(4), 0; 0, t201 * t211 + t203 * t210, 0, t206, 0; 0, -t207, 0, 0, 0; 0, -t208, 0, 0, 0; 0, -t212, 0, 0, 0; 0, t206 * t200, 0, t201 * t207 + (-t199 * t201 + t200 * t213) * qJD(4), 0; 0, t206 * t199, 0, t201 * t208 + (t199 * t213 + t200 * t201) * qJD(4), 0; 0, t201 * t210 - t203 * t211, 0, t205, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end