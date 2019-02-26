% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRPR7
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:14
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRRPR7_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR7_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR7_jacobig_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:14:05
% EndTime: 2019-02-26 20:14:05
% DurationCPUTime: 0.06s
% Computational Cost: add. (16->14), mult. (48->33), div. (0->0), fcn. (72->10), ass. (0->18)
t154 = sin(pkin(12));
t156 = sin(pkin(6));
t168 = t154 * t156;
t155 = sin(pkin(7));
t167 = t155 * t156;
t157 = cos(pkin(12));
t166 = t157 * t156;
t159 = cos(pkin(6));
t161 = sin(qJ(2));
t165 = t159 * t161;
t163 = cos(qJ(2));
t164 = t159 * t163;
t162 = cos(qJ(3));
t160 = sin(qJ(3));
t158 = cos(pkin(7));
t153 = -t154 * t164 - t157 * t161;
t152 = -t154 * t161 + t157 * t164;
t1 = [0, t168, -t153 * t155 + t158 * t168 (-t154 * t165 + t157 * t163) * t160 + (-t153 * t158 - t154 * t167) * t162, 0, 0; 0, -t166, -t152 * t155 - t158 * t166 (t154 * t163 + t157 * t165) * t160 + (-t152 * t158 + t155 * t166) * t162, 0, 0; 0, t159, t159 * t158 - t163 * t167, -t159 * t155 * t162 + (-t158 * t162 * t163 + t160 * t161) * t156, 0, 0;];
Jg_rot  = t1;
