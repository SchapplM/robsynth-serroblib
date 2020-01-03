% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPRRP8
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RPRRP8_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP8_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP8_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPRRP8_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRP8_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP8_jacobiaD_transl_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:47:46
	% EndTime: 2019-12-31 18:47:46
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:47:46
	% EndTime: 2019-12-31 18:47:46
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:47:46
	% EndTime: 2019-12-31 18:47:46
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->6), mult. (20->10), div. (0->0), fcn. (12->2), ass. (0->5)
	t10 = -pkin(1) - r_i_i_C(1);
	t9 = r_i_i_C(3) + qJ(2);
	t8 = cos(qJ(1));
	t7 = sin(qJ(1));
	t1 = [t8 * qJD(2) + (t10 * t8 - t9 * t7) * qJD(1), qJD(1) * t8, 0, 0, 0; t7 * qJD(2) + (t10 * t7 + t9 * t8) * qJD(1), qJD(1) * t7, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:47:46
	% EndTime: 2019-12-31 18:47:46
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (36->15), mult. (88->23), div. (0->0), fcn. (74->4), ass. (0->14)
	t87 = -pkin(1) - pkin(2);
	t78 = sin(qJ(1));
	t86 = qJD(1) * t78;
	t80 = cos(qJ(1));
	t85 = qJD(1) * t80;
	t84 = qJD(3) * t78;
	t83 = qJD(3) * t80;
	t77 = sin(qJ(3));
	t79 = cos(qJ(3));
	t71 = -t77 * t84 - t79 * t83 + (t77 * t78 + t79 * t80) * qJD(1);
	t72 = (t84 - t86) * t79 + (-t83 + t85) * t77;
	t82 = -t72 * r_i_i_C(1) - t71 * r_i_i_C(2);
	t81 = t71 * r_i_i_C(1) - t72 * r_i_i_C(2);
	t1 = [t80 * qJD(2) + (-qJ(2) * t78 + t87 * t80) * qJD(1) - t81, t85, t81, 0, 0; t78 * qJD(2) + (qJ(2) * t80 + t87 * t78) * qJD(1) - t82, t86, t82, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:47:46
	% EndTime: 2019-12-31 18:47:46
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (123->23), mult. (308->39), div. (0->0), fcn. (276->6), ass. (0->21)
	t80 = qJD(1) - qJD(3);
	t68 = sin(qJ(3));
	t69 = sin(qJ(1));
	t70 = cos(qJ(3));
	t71 = cos(qJ(1));
	t47 = -t69 * t68 - t71 * t70;
	t45 = t80 * t47;
	t48 = t71 * t68 - t69 * t70;
	t46 = t80 * t48;
	t56 = sin(qJ(4));
	t57 = cos(qJ(4));
	t60 = r_i_i_C(1) * t57 - r_i_i_C(2) * t56 + pkin(3);
	t72 = pkin(7) + r_i_i_C(3);
	t77 = (r_i_i_C(1) * t56 + r_i_i_C(2) * t57) * qJD(4);
	t79 = -t72 * t45 - t60 * t46 - t47 * t77;
	t78 = -t60 * t45 + t72 * t46 + t48 * t77;
	t76 = qJ(2) * qJD(1);
	t73 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
	t67 = qJD(4) * t56;
	t66 = qJD(4) * t57;
	t1 = [-t69 * t76 + t73 * t71 - t78, qJD(1) * t71, t78, (-t46 * t57 - t47 * t67) * r_i_i_C(2) + (-t46 * t56 + t47 * t66) * r_i_i_C(1), 0; t73 * t69 + t71 * t76 - t79, qJD(1) * t69, t79, (t45 * t57 - t48 * t67) * r_i_i_C(2) + (t45 * t56 + t48 * t66) * r_i_i_C(1), 0; 0, 0, 0, t77, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:47:47
	% EndTime: 2019-12-31 18:47:47
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (221->26), mult. (544->37), div. (0->0), fcn. (506->6), ass. (0->26)
	t275 = sin(qJ(3));
	t276 = sin(qJ(1));
	t277 = cos(qJ(3));
	t278 = cos(qJ(1));
	t250 = -t276 * t275 - t278 * t277;
	t287 = qJD(1) - qJD(3);
	t248 = t287 * t250;
	t251 = t278 * t275 - t276 * t277;
	t249 = t287 * t251;
	t259 = sin(qJ(4));
	t260 = cos(qJ(4));
	t274 = r_i_i_C(3) + qJ(5);
	t280 = pkin(4) + r_i_i_C(1);
	t281 = t280 * t259 - t274 * t260;
	t262 = t281 * qJD(4) - t259 * qJD(5);
	t267 = t274 * t259 + t280 * t260;
	t264 = pkin(3) + t267;
	t279 = pkin(7) + r_i_i_C(2);
	t291 = -t279 * t248 - t264 * t249 - t262 * t250;
	t290 = -t264 * t248 + t279 * t249 + t262 * t251;
	t284 = qJD(1) * t276;
	t283 = qJD(1) * t278;
	t282 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
	t273 = qJD(4) * t260;
	t263 = t267 * qJD(4) - qJD(5) * t260;
	t1 = [-qJ(2) * t284 + t282 * t278 - t290, t283, t290, -t249 * t281 + t263 * t250, t249 * t259 - t250 * t273; qJ(2) * t283 + t282 * t276 - t291, t284, t291, t248 * t281 + t263 * t251, -t248 * t259 - t251 * t273; 0, 0, 0, t262, -qJD(4) * t259;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end