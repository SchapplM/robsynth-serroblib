% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPPPRR5
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:32
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPPRR5_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR5_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPPRR5_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPPRR5_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:32:36
	% EndTime: 2019-10-09 23:32:36
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:32:36
	% EndTime: 2019-10-09 23:32:36
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:32:36
	% EndTime: 2019-10-09 23:32:36
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (8->6), mult. (20->10), div. (0->0), fcn. (12->2), ass. (0->5)
	t10 = -pkin(1) + r_i_i_C(2);
	t9 = r_i_i_C(3) + qJ(2);
	t8 = cos(qJ(1));
	t7 = sin(qJ(1));
	t1 = [t8 * qJD(2) + (t10 * t8 - t9 * t7) * qJD(1), qJD(1) * t8, 0, 0, 0, 0; t7 * qJD(2) + (t10 * t7 + t9 * t8) * qJD(1), qJD(1) * t7, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:32:36
	% EndTime: 2019-10-09 23:32:36
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (13->9), mult. (28->12), div. (0->0), fcn. (18->2), ass. (0->7)
	t12 = r_i_i_C(1) + qJ(2);
	t8 = sin(qJ(1));
	t11 = qJD(1) * t8;
	t10 = -pkin(1) - r_i_i_C(3) - qJ(3);
	t9 = cos(qJ(1));
	t7 = qJD(1) * t9;
	t1 = [t9 * qJD(2) - t8 * qJD(3) + (t10 * t9 - t12 * t8) * qJD(1), t7, -t11, 0, 0, 0; t8 * qJD(2) + t9 * qJD(3) + (t10 * t8 + t12 * t9) * qJD(1), t11, t7, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:32:36
	% EndTime: 2019-10-09 23:32:36
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (19->15), mult. (44->22), div. (0->0), fcn. (32->4), ass. (0->11)
	t55 = -pkin(1) - qJ(3);
	t54 = pkin(3) + qJ(2);
	t51 = sin(qJ(1));
	t53 = qJD(1) * t51;
	t52 = cos(qJ(1));
	t50 = cos(pkin(9));
	t49 = sin(pkin(9));
	t48 = qJD(1) * t52;
	t47 = (-t49 * t51 + t50 * t52) * qJD(1);
	t46 = (-t49 * t52 - t50 * t51) * qJD(1);
	t1 = [t46 * r_i_i_C(1) - t47 * r_i_i_C(2) + t52 * qJD(2) - t51 * qJD(3) + (-t54 * t51 + t55 * t52) * qJD(1), t48, -t53, 0, 0, 0; t47 * r_i_i_C(1) + t46 * r_i_i_C(2) + t51 * qJD(2) + t52 * qJD(3) + (t55 * t51 + t54 * t52) * qJD(1), t53, t48, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:32:36
	% EndTime: 2019-10-09 23:32:36
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (54->27), mult. (146->43), div. (0->0), fcn. (122->6), ass. (0->20)
	t56 = pkin(7) + r_i_i_C(3);
	t55 = -pkin(1) - qJ(3);
	t54 = pkin(3) + qJ(2);
	t46 = sin(qJ(1));
	t53 = qJD(1) * t46;
	t45 = sin(qJ(5));
	t52 = qJD(5) * t45;
	t47 = cos(qJ(5));
	t51 = qJD(5) * t47;
	t43 = sin(pkin(9));
	t44 = cos(pkin(9));
	t48 = cos(qJ(1));
	t40 = t43 * t48 + t46 * t44;
	t39 = -t46 * t43 + t44 * t48;
	t50 = r_i_i_C(1) * t47 - r_i_i_C(2) * t45 + pkin(4);
	t49 = (-r_i_i_C(1) * t45 - r_i_i_C(2) * t47) * qJD(5);
	t42 = qJD(1) * t48;
	t38 = t39 * qJD(1);
	t37 = t40 * qJD(1);
	t1 = [t48 * qJD(2) - t46 * qJD(3) + t56 * t38 + t39 * t49 - t50 * t37 + (-t54 * t46 + t55 * t48) * qJD(1), t42, -t53, 0, (-t38 * t47 + t40 * t52) * r_i_i_C(2) + (-t38 * t45 - t40 * t51) * r_i_i_C(1), 0; t46 * qJD(2) + t48 * qJD(3) + t56 * t37 + t40 * t49 + t50 * t38 + (t55 * t46 + t54 * t48) * qJD(1), t53, t42, 0, (-t37 * t47 - t39 * t52) * r_i_i_C(2) + (-t37 * t45 + t39 * t51) * r_i_i_C(1), 0; 0, 0, 0, 0, t49, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:32:37
	% EndTime: 2019-10-09 23:32:38
	% DurationCPUTime: 0.36s
	% Computational Cost: add. (176->54), mult. (490->85), div. (0->0), fcn. (467->8), ass. (0->40)
	t224 = cos(qJ(5));
	t221 = sin(qJ(6));
	t223 = cos(qJ(6));
	t232 = t223 * r_i_i_C(1) - t221 * r_i_i_C(2) + pkin(5);
	t222 = sin(qJ(5));
	t248 = pkin(8) + r_i_i_C(3);
	t239 = t248 * t222;
	t252 = t232 * t224 + t239;
	t253 = t252 * qJD(5);
	t251 = -pkin(1) - qJ(3);
	t250 = pkin(3) + qJ(2);
	t247 = sin(qJ(1));
	t246 = sin(pkin(9));
	t225 = cos(qJ(1));
	t219 = qJD(1) * t225;
	t245 = qJD(5) * t222;
	t244 = qJD(5) * t224;
	t220 = cos(pkin(9));
	t216 = t247 * t220 + t225 * t246;
	t243 = qJD(6) * t216;
	t242 = qJD(6) * t222;
	t241 = qJD(6) * t223;
	t240 = qJD(6) * t224;
	t238 = t248 * t224;
	t237 = qJD(1) * t247;
	t236 = t247 * t246;
	t235 = t221 * r_i_i_C(1) + t223 * r_i_i_C(2);
	t213 = t216 * qJD(1);
	t234 = -t216 * t240 + t213;
	t214 = -qJD(1) * t236 + t220 * t219;
	t215 = t225 * t220 - t236;
	t233 = t215 * t240 - t214;
	t231 = qJD(6) * t235;
	t230 = -t213 * t224 - t215 * t245 + t243;
	t229 = qJD(6) * t215 - t214 * t224 + t216 * t245;
	t227 = -t232 * t222 + t238;
	t226 = t227 * qJD(5) - t224 * t231;
	t212 = t234 * t221 - t229 * t223;
	t211 = t229 * t221 + t234 * t223;
	t1 = [(t214 * t221 + t216 * t241) * r_i_i_C(1) + (t214 * t223 - t221 * t243) * r_i_i_C(2) + t214 * pkin(7) - t247 * qJD(3) + t225 * qJD(2) + (-pkin(4) - t252) * t213 + (t251 * t225 - t250 * t247) * qJD(1) + t226 * t215, t219, -t237, 0, t227 * t214 + (t235 * t242 - t253) * t216, t211 * r_i_i_C(1) - t212 * r_i_i_C(2); t247 * qJD(2) + t213 * pkin(7) + t212 * r_i_i_C(1) + t211 * r_i_i_C(2) + t225 * qJD(3) + (-pkin(5) * t222 + t238) * t216 * qJD(5) + (pkin(5) * t224 + pkin(4) + t239) * t214 + (t250 * t225 + t251 * t247) * qJD(1), t237, t219, 0, t227 * t213 + (-t222 * t231 + t253) * t215, (t233 * r_i_i_C(1) + t230 * r_i_i_C(2)) * t223 + (t230 * r_i_i_C(1) - t233 * r_i_i_C(2)) * t221; 0, 0, 0, 0, t226, (t221 * t242 - t223 * t244) * r_i_i_C(2) + (-t221 * t244 - t222 * t241) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end