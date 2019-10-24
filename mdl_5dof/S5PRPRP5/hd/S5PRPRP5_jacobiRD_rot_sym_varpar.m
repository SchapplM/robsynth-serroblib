% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PRPRP5
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:24
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5PRPRP5_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP5_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP5_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPRP5_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP5_jacobiRD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:24:49
	% EndTime: 2019-10-24 10:24:49
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:24:49
	% EndTime: 2019-10-24 10:24:49
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:24:48
	% EndTime: 2019-10-24 10:24:49
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
	% StartTime: 2019-10-24 10:24:49
	% EndTime: 2019-10-24 10:24:49
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (5->5), mult. (21->12), div. (0->0), fcn. (21->6), ass. (0->9)
	t133 = qJD(2) * sin(qJ(2));
	t132 = qJD(2) * cos(qJ(2));
	t125 = sin(pkin(7));
	t131 = t125 * t132;
	t127 = cos(pkin(7));
	t130 = t127 * t132;
	t126 = cos(pkin(8));
	t124 = sin(pkin(8));
	t1 = [0, -t126 * t130, 0, 0, 0; 0, -t126 * t131, 0, 0, 0; 0, -t126 * t133, 0, 0, 0; 0, t124 * t130, 0, 0, 0; 0, t124 * t131, 0, 0, 0; 0, t124 * t133, 0, 0, 0; 0, -t127 * t133, 0, 0, 0; 0, -t125 * t133, 0, 0, 0; 0, t132, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:24:49
	% EndTime: 2019-10-24 10:24:49
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (46->16), mult. (77->36), div. (0->0), fcn. (77->6), ass. (0->18)
	t161 = sin(pkin(7));
	t164 = cos(qJ(2));
	t174 = t161 * t164;
	t162 = cos(pkin(7));
	t173 = t162 * t164;
	t163 = sin(qJ(2));
	t172 = qJD(2) * t163;
	t171 = qJD(2) * t164;
	t170 = qJD(4) * t163;
	t169 = qJD(4) * t164;
	t168 = t161 * t172;
	t167 = t162 * t172;
	t160 = pkin(8) + qJ(4);
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
	% StartTime: 2019-10-24 10:24:49
	% EndTime: 2019-10-24 10:24:50
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (46->17), mult. (77->36), div. (0->0), fcn. (77->6), ass. (0->18)
	t216 = sin(pkin(7));
	t219 = cos(qJ(2));
	t229 = t216 * t219;
	t217 = cos(pkin(7));
	t228 = t217 * t219;
	t218 = sin(qJ(2));
	t227 = qJD(2) * t218;
	t226 = qJD(2) * t219;
	t225 = qJD(4) * t218;
	t224 = qJD(4) * t219;
	t223 = t216 * t227;
	t222 = t217 * t227;
	t215 = pkin(8) + qJ(4);
	t213 = sin(t215);
	t214 = cos(t215);
	t221 = -t213 * t225 + t214 * t226;
	t220 = -t213 * t226 - t214 * t225;
	t1 = [0, -t221 * t217, 0, t213 * t222 + (-t213 * t216 - t214 * t228) * qJD(4), 0; 0, -t221 * t216, 0, t213 * t223 + (t213 * t217 - t214 * t229) * qJD(4), 0; 0, -t213 * t224 - t214 * t227, 0, t220, 0; 0, -t222, 0, 0, 0; 0, -t223, 0, 0, 0; 0, t226, 0, 0, 0; 0, t220 * t217, 0, -t214 * t222 + (-t213 * t228 + t214 * t216) * qJD(4), 0; 0, t220 * t216, 0, -t214 * t223 + (-t213 * t229 - t214 * t217) * qJD(4), 0; 0, -t213 * t227 + t214 * t224, 0, t221, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end