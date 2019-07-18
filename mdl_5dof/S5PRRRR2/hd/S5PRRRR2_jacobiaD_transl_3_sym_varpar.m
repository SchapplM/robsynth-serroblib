% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S5PRRRR2
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,d5]';
%
% Output:
% JaD_transl [3x5]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:30
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5PRRRR2_jacobiaD_transl_3_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_jacobiaD_transl_3_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR2_jacobiaD_transl_3_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRRRR2_jacobiaD_transl_3_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR2_jacobiaD_transl_3_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:30:34
% EndTime: 2019-07-18 13:30:34
% DurationCPUTime: 0.04s
% Computational Cost: add. (22->6), mult. (20->9), div. (0->0), fcn. (10->4), ass. (0->8)
t43 = pkin(2) * qJD(2);
t40 = qJ(2) + qJ(3);
t37 = sin(t40);
t38 = cos(t40);
t39 = qJD(2) + qJD(3);
t42 = (-r_i_i_C(1) * t38 + r_i_i_C(2) * t37) * t39;
t41 = (-r_i_i_C(1) * t37 - r_i_i_C(2) * t38) * t39;
t1 = [0, -cos(qJ(2)) * t43 + t42, t42, 0, 0; 0, -sin(qJ(2)) * t43 + t41, t41, 0, 0; 0, 0, 0, 0, 0;];
JaD_transl  = t1;
