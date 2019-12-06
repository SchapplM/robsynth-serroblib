% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PRPRR7
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5PRPRR7_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR7_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR7_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPRR7_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR7_jacobiRD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:01:10
	% EndTime: 2019-12-05 16:01:11
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:01:10
	% EndTime: 2019-12-05 16:01:10
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:01:11
	% EndTime: 2019-12-05 16:01:11
	% DurationCPUTime: 0.05s
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
	% StartTime: 2019-12-05 16:01:11
	% EndTime: 2019-12-05 16:01:11
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (2->2), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->5)
	t103 = qJD(2) * sin(qJ(2));
	t102 = qJD(2) * cos(qJ(2));
	t99 = cos(pkin(8));
	t98 = sin(pkin(8));
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, t99 * t102, 0, 0, 0; 0, t98 * t102, 0, 0, 0; 0, t103, 0, 0, 0; 0, -t99 * t103, 0, 0, 0; 0, -t98 * t103, 0, 0, 0; 0, t102, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:01:12
	% EndTime: 2019-12-05 16:01:12
	% DurationCPUTime: 0.05s
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
	t146 = sin(pkin(8));
	t155 = t146 * t158;
	t147 = cos(pkin(8));
	t154 = t147 * t158;
	t153 = t148 * t156 + t150 * t159;
	t152 = -t148 * t159 + t150 * t156;
	t1 = [0, t152 * t147, 0, t150 * t154 + (-t146 * t150 - t147 * t161) * qJD(4), 0; 0, t152 * t146, 0, t150 * t155 + (-t146 * t161 + t147 * t150) * qJD(4), 0; 0, t148 * t158 + t150 * t157, 0, t153, 0; 0, -t153 * t147, 0, -t148 * t154 + (t146 * t148 - t147 * t160) * qJD(4), 0; 0, -t153 * t146, 0, -t148 * t155 + (-t146 * t160 - t147 * t148) * qJD(4), 0; 0, -t148 * t157 + t150 * t158, 0, t152, 0; 0, -t154, 0, 0, 0; 0, -t155, 0, 0, 0; 0, -t159, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:01:12
	% EndTime: 2019-12-05 16:01:12
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (99->17), mult. (117->28), div. (0->0), fcn. (117->6), ass. (0->25)
	t209 = qJD(4) + qJD(5);
	t213 = sin(qJ(2));
	t224 = t209 * t213;
	t214 = cos(qJ(2));
	t223 = t209 * t214;
	t222 = qJD(2) * t213;
	t221 = qJD(2) * t214;
	t210 = qJ(4) + qJ(5);
	t207 = sin(t210);
	t220 = t207 * t224;
	t208 = cos(t210);
	t219 = t208 * t224;
	t211 = sin(pkin(8));
	t218 = t211 * t221;
	t212 = cos(pkin(8));
	t217 = t212 * t221;
	t216 = -t209 * t211 + t217;
	t215 = t209 * t212 + t218;
	t206 = t207 * t223 + t208 * t222;
	t205 = -t207 * t222 + t208 * t223;
	t204 = t215 * t208 - t211 * t220;
	t203 = -t215 * t207 - t211 * t219;
	t202 = t216 * t208 - t212 * t220;
	t201 = -t216 * t207 - t212 * t219;
	t1 = [0, t205 * t212, 0, t202, t202; 0, t205 * t211, 0, t204, t204; 0, t207 * t221 + t219, 0, t206, t206; 0, -t206 * t212, 0, t201, t201; 0, -t206 * t211, 0, t203, t203; 0, t208 * t221 - t220, 0, t205, t205; 0, -t217, 0, 0, 0; 0, -t218, 0, 0, 0; 0, -t222, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end