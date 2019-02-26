% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PPRRPR2
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PPRRPR2_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR2_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_jacobigD_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:40:57
% EndTime: 2019-02-26 19:40:57
% DurationCPUTime: 0.08s
% Computational Cost: add. (10->10), mult. (40->28), div. (0->0), fcn. (44->10), ass. (0->14)
t156 = sin(pkin(11));
t162 = cos(pkin(6));
t167 = t156 * t162;
t157 = sin(pkin(7));
t158 = sin(pkin(6));
t166 = t158 * t157;
t160 = cos(pkin(11));
t165 = t160 * t162;
t164 = cos(qJ(3));
t163 = sin(qJ(3));
t161 = cos(pkin(7));
t159 = cos(pkin(12));
t155 = sin(pkin(12));
t1 = [0, 0, 0 ((-t155 * t167 + t159 * t160) * t164 + ((-t155 * t160 - t159 * t167) * t161 + t156 * t166) * t163) * qJD(3), 0, 0; 0, 0, 0 ((t155 * t165 + t156 * t159) * t164 + ((-t155 * t156 + t159 * t165) * t161 - t160 * t166) * t163) * qJD(3), 0, 0; 0, 0, 0 (t157 * t162 * t163 + (t159 * t161 * t163 + t155 * t164) * t158) * qJD(3), 0, 0;];
JgD_rot  = t1;
