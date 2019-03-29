% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S5RRRRR2
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
%
% Output:
% JaD_transl [3x5]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-29 15:26
% Revision: 932832b1be1be80f59b7f1a581a1a8f328bdb39d (2019-03-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRRRR2_jacobiaD_transl_3_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(3,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_jacobiaD_transl_3_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR2_jacobiaD_transl_3_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRRR2_jacobiaD_transl_3_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_jacobiaD_transl_3_sym_varpar: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-29 15:26:51
% EndTime: 2019-03-29 15:26:51
% DurationCPUTime: 0.08s
% Computational Cost: add. (69->16), mult. (88->32), div. (0->0), fcn. (56->6), ass. (0->18)
t36 = qJ(1) + qJ(2);
t33 = sin(t36);
t34 = cos(t36);
t38 = cos(qJ(3));
t47 = qJD(3) * t38;
t35 = qJD(1) + qJD(2);
t37 = sin(qJ(3));
t51 = t35 * t37;
t53 = t33 * t51 - t34 * t47;
t52 = t33 * t47 + t34 * t51;
t50 = t35 * t38;
t49 = pkin(1) * qJD(1);
t48 = qJD(3) * t37;
t44 = t33 * t48;
t41 = t33 * t50 + t34 * t48;
t40 = r_i_i_C(1) * t44 + (-r_i_i_C(1) * t34 * t38 - r_i_i_C(3) * t33) * t35 + t52 * r_i_i_C(2);
t39 = t35 * t34 * r_i_i_C(3) - t41 * r_i_i_C(1) + t53 * r_i_i_C(2);
t1 = [-cos(qJ(1)) * t49 + t40, t40, t53 * r_i_i_C(1) + t41 * r_i_i_C(2), 0, 0; -sin(qJ(1)) * t49 + t39, t39 (-t34 * t50 + t44) * r_i_i_C(2) - t52 * r_i_i_C(1), 0, 0; 0, 0 (-r_i_i_C(1) * t37 - r_i_i_C(2) * t38) * qJD(3), 0, 0;];
JaD_transl  = t1;
