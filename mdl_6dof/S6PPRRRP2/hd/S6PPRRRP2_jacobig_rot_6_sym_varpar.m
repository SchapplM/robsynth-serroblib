% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRRRP2
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:42
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PPRRRP2_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP2_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP2_jacobig_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:42:05
% EndTime: 2019-02-26 19:42:05
% DurationCPUTime: 0.06s
% Computational Cost: add. (33->20), mult. (98->46), div. (0->0), fcn. (140->12), ass. (0->28)
t162 = sin(pkin(11));
t168 = cos(pkin(6));
t180 = t162 * t168;
t163 = sin(pkin(7));
t164 = sin(pkin(6));
t179 = t163 * t164;
t178 = t163 * t168;
t167 = cos(pkin(7));
t177 = t164 * t167;
t165 = cos(pkin(12));
t176 = t165 * t167;
t166 = cos(pkin(11));
t175 = t166 * t168;
t161 = sin(pkin(12));
t157 = -t162 * t161 + t165 * t175;
t174 = -t157 * t167 + t166 * t179;
t159 = -t166 * t161 - t165 * t180;
t173 = t159 * t167 + t162 * t179;
t172 = cos(qJ(3));
t171 = cos(qJ(4));
t170 = sin(qJ(3));
t169 = sin(qJ(4));
t160 = -t161 * t180 + t166 * t165;
t158 = t161 * t175 + t162 * t165;
t156 = -t165 * t179 + t168 * t167;
t155 = -t159 * t163 + t162 * t177;
t154 = -t157 * t163 - t166 * t177;
t1 = [0, 0, t155, t160 * t170 - t173 * t172 (t160 * t172 + t173 * t170) * t169 - t155 * t171, 0; 0, 0, t154, t158 * t170 + t174 * t172 (t158 * t172 - t174 * t170) * t169 - t154 * t171, 0; 0, 0, t156, -t172 * t178 + (t161 * t170 - t172 * t176) * t164 (t170 * t178 + (t161 * t172 + t170 * t176) * t164) * t169 - t156 * t171, 0;];
Jg_rot  = t1;
