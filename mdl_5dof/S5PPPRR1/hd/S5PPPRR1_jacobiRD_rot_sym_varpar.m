% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PPPRR1
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
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 14:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5PPPRR1_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR1_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR1_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPPRR1_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR1_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 14:58:24
	% EndTime: 2019-12-05 14:58:24
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 14:58:24
	% EndTime: 2019-12-05 14:58:24
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 14:58:24
	% EndTime: 2019-12-05 14:58:24
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 14:58:24
	% EndTime: 2019-12-05 14:58:24
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 14:58:24
	% EndTime: 2019-12-05 14:58:24
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (15->6), mult. (24->17), div. (0->0), fcn. (24->6), ass. (0->10)
	t70 = sin(pkin(7));
	t71 = cos(pkin(8));
	t75 = t70 * t71;
	t72 = cos(pkin(7));
	t74 = t71 * t72;
	t73 = qJD(4) * sin(pkin(8));
	t68 = pkin(9) + qJ(4);
	t67 = cos(t68);
	t66 = sin(t68);
	t1 = [0, 0, 0, (-t66 * t70 - t67 * t74) * qJD(4), 0; 0, 0, 0, (t66 * t72 - t67 * t75) * qJD(4), 0; 0, 0, 0, -t67 * t73, 0; 0, 0, 0, (t66 * t74 - t67 * t70) * qJD(4), 0; 0, 0, 0, (t66 * t75 + t67 * t72) * qJD(4), 0; 0, 0, 0, t66 * t73, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 14:58:24
	% EndTime: 2019-12-05 14:58:25
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (82->24), mult. (140->57), div. (0->0), fcn. (148->8), ass. (0->27)
	t219 = sin(pkin(8));
	t223 = sin(qJ(5));
	t233 = t219 * t223;
	t224 = cos(qJ(5));
	t232 = t219 * t224;
	t220 = sin(pkin(7));
	t221 = cos(pkin(8));
	t231 = t220 * t221;
	t218 = pkin(9) + qJ(4);
	t216 = sin(t218);
	t222 = cos(pkin(7));
	t230 = t222 * t216;
	t217 = cos(t218);
	t229 = t222 * t217;
	t228 = qJD(4) * t217;
	t227 = qJD(5) * t223;
	t226 = qJD(5) * t224;
	t225 = t219 * qJD(4) * t216;
	t215 = t220 * t216 + t221 * t229;
	t213 = t217 * t231 - t230;
	t214 = t220 * t217 - t221 * t230;
	t212 = -t216 * t231 - t229;
	t211 = t215 * qJD(4);
	t210 = t214 * qJD(4);
	t209 = t213 * qJD(4);
	t208 = t212 * qJD(4);
	t1 = [0, 0, 0, -t211 * t224 - t214 * t227, -t210 * t223 + (-t215 * t224 - t222 * t233) * qJD(5); 0, 0, 0, -t209 * t224 - t212 * t227, -t208 * t223 + (-t213 * t224 - t220 * t233) * qJD(5); 0, 0, 0, (t216 * t227 - t224 * t228) * t219, t223 * t225 + (-t217 * t232 + t221 * t223) * qJD(5); 0, 0, 0, t211 * t223 - t214 * t226, -t210 * t224 + (t215 * t223 - t222 * t232) * qJD(5); 0, 0, 0, t209 * t223 - t212 * t226, -t208 * t224 + (t213 * t223 - t220 * t232) * qJD(5); 0, 0, 0, (t216 * t226 + t223 * t228) * t219, t224 * t225 + (t217 * t233 + t221 * t224) * qJD(5); 0, 0, 0, t210, 0; 0, 0, 0, t208, 0; 0, 0, 0, -t225, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end