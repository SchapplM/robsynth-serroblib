% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPPRR3
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPPPRR3_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_jacobiR_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:23:51
% EndTime: 2019-02-26 20:23:51
% DurationCPUTime: 0.02s
% Computational Cost: add. (27->6), mult. (28->8), div. (0->0), fcn. (50->6), ass. (0->14)
t39 = cos(qJ(1));
t38 = sin(qJ(1));
t32 = sin(pkin(9));
t33 = cos(pkin(9));
t24 = -t38 * t32 - t39 * t33;
t31 = pkin(10) + qJ(5);
t29 = sin(t31);
t37 = t24 * t29;
t30 = cos(t31);
t36 = t24 * t30;
t25 = t39 * t32 - t38 * t33;
t35 = t25 * t29;
t34 = t25 * t30;
t1 = [t34, 0, 0, 0, t37, 0; -t36, 0, 0, 0, t35, 0; 0, 0, 0, 0, -t30, 0; -t35, 0, 0, 0, t36, 0; t37, 0, 0, 0, t34, 0; 0, 0, 0, 0, t29, 0; t24, 0, 0, 0, 0, 0; t25, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
JR_rot  = t1;
