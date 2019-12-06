% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PRPRR5
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
% Datum: 2019-12-05 15:55
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5PRPRR5_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR5_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR5_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPRR5_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR5_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:55:21
	% EndTime: 2019-12-05 15:55:21
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:55:21
	% EndTime: 2019-12-05 15:55:21
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:55:21
	% EndTime: 2019-12-05 15:55:21
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
	% StartTime: 2019-12-05 15:55:21
	% EndTime: 2019-12-05 15:55:22
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (5->5), mult. (21->12), div. (0->0), fcn. (21->6), ass. (0->9)
	t133 = qJD(2) * sin(qJ(2));
	t132 = qJD(2) * cos(qJ(2));
	t125 = sin(pkin(8));
	t131 = t125 * t132;
	t127 = cos(pkin(8));
	t130 = t127 * t132;
	t126 = cos(pkin(9));
	t124 = sin(pkin(9));
	t1 = [0, -t126 * t130, 0, 0, 0; 0, -t126 * t131, 0, 0, 0; 0, -t126 * t133, 0, 0, 0; 0, t124 * t130, 0, 0, 0; 0, t124 * t131, 0, 0, 0; 0, t124 * t133, 0, 0, 0; 0, -t127 * t133, 0, 0, 0; 0, -t125 * t133, 0, 0, 0; 0, t132, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:55:22
	% EndTime: 2019-12-05 15:55:22
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (46->16), mult. (77->36), div. (0->0), fcn. (77->6), ass. (0->18)
	t161 = sin(pkin(8));
	t164 = cos(qJ(2));
	t174 = t161 * t164;
	t162 = cos(pkin(8));
	t173 = t162 * t164;
	t163 = sin(qJ(2));
	t172 = qJD(2) * t163;
	t171 = qJD(2) * t164;
	t170 = qJD(4) * t163;
	t169 = qJD(4) * t164;
	t168 = t161 * t172;
	t167 = t162 * t172;
	t160 = pkin(9) + qJ(4);
	t158 = sin(t160);
	t159 = cos(t160);
	t166 = t158 * t170 - t159 * t171;
	t165 = t158 * t171 + t159 * t170;
	t1 = [0, t166 * t162, 0, t158 * t167 + (-t158 * t161 - t159 * t173) * qJD(4), 0; 0, t166 * t161, 0, t158 * t168 + (t158 * t162 - t159 * t174) * qJD(4), 0; 0, -t158 * t169 - t159 * t172, 0, -t165, 0; 0, t165 * t162, 0, t159 * t167 + (t158 * t173 - t159 * t161) * qJD(4), 0; 0, t165 * t161, 0, t159 * t168 + (t158 * t174 + t159 * t162) * qJD(4), 0; 0, t158 * t172 - t159 * t169, 0, t166, 0; 0, -t167, 0, 0, 0; 0, -t168, 0, 0, 0; 0, t171, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:55:22
	% EndTime: 2019-12-05 15:55:22
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (142->16), mult. (117->28), div. (0->0), fcn. (117->6), ass. (0->25)
	t214 = qJD(4) + qJD(5);
	t217 = sin(qJ(2));
	t229 = t214 * t217;
	t218 = cos(qJ(2));
	t228 = t214 * t218;
	t227 = qJD(2) * t217;
	t226 = qJD(2) * t218;
	t213 = pkin(9) + qJ(4) + qJ(5);
	t211 = sin(t213);
	t225 = t211 * t228;
	t212 = cos(t213);
	t224 = t212 * t228;
	t215 = sin(pkin(8));
	t223 = t215 * t227;
	t216 = cos(pkin(8));
	t222 = t216 * t227;
	t221 = -t214 * t215 + t222;
	t220 = t214 * t216 + t223;
	t210 = t211 * t229 - t212 * t226;
	t219 = t211 * t226 + t212 * t229;
	t208 = t221 * t212 + t216 * t225;
	t207 = t221 * t211 - t216 * t224;
	t206 = t220 * t212 + t215 * t225;
	t205 = t220 * t211 - t215 * t224;
	t1 = [0, t210 * t216, 0, t207, t207; 0, t210 * t215, 0, t205, t205; 0, -t212 * t227 - t225, 0, -t219, -t219; 0, t219 * t216, 0, t208, t208; 0, t219 * t215, 0, t206, t206; 0, t211 * t227 - t224, 0, t210, t210; 0, -t222, 0, 0, 0; 0, -t223, 0, 0, 0; 0, t226, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end