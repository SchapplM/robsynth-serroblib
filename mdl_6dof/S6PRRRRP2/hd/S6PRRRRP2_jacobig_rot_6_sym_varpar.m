% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRP2
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:15
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRRRP2_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP2_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP2_jacobig_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:15:43
% EndTime: 2019-02-26 20:15:43
% DurationCPUTime: 0.03s
% Computational Cost: add. (18->11), mult. (31->20), div. (0->0), fcn. (52->8), ass. (0->17)
t156 = sin(pkin(11));
t157 = sin(pkin(6));
t166 = t156 * t157;
t161 = cos(qJ(2));
t165 = t157 * t161;
t158 = cos(pkin(11));
t164 = t158 * t157;
t159 = cos(pkin(6));
t160 = sin(qJ(2));
t163 = t159 * t160;
t162 = t159 * t161;
t155 = qJ(3) + qJ(4);
t154 = cos(t155);
t153 = sin(t155);
t152 = t156 * t162 + t158 * t160;
t151 = t156 * t160 - t158 * t162;
t1 = [0, t166, t152, t152 (-t156 * t163 + t158 * t161) * t153 - t154 * t166, 0; 0, -t164, t151, t151 (t156 * t161 + t158 * t163) * t153 + t154 * t164, 0; 0, t159, -t165, -t165, t157 * t160 * t153 - t159 * t154, 0;];
Jg_rot  = t1;
