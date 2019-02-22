% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPRPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 10:10
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPPRPR4_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR4_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR4_jacobiR_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:10:27
% EndTime: 2019-02-22 10:10:27
% DurationCPUTime: 0.02s
% Computational Cost: add. (27->6), mult. (28->8), div. (0->0), fcn. (50->6), ass. (0->14)
t43 = cos(qJ(1));
t42 = sin(qJ(1));
t36 = sin(pkin(9));
t37 = cos(pkin(9));
t28 = -t42 * t36 - t43 * t37;
t35 = qJ(4) + pkin(10);
t33 = sin(t35);
t41 = t28 * t33;
t34 = cos(t35);
t40 = t28 * t34;
t29 = t43 * t36 - t42 * t37;
t39 = t29 * t33;
t38 = t29 * t34;
t1 = [t38, 0, 0, t41, 0, 0; -t40, 0, 0, t39, 0, 0; 0, 0, 0, -t34, 0, 0; -t39, 0, 0, t40, 0, 0; t41, 0, 0, t38, 0, 0; 0, 0, 0, t33, 0, 0; t28, 0, 0, 0, 0, 0; t29, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
JR_rot  = t1;
