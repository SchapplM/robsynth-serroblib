% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRRR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:02
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRRR7_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR7_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR7_jacobiR_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:02:00
% EndTime: 2019-02-22 11:02:00
% DurationCPUTime: 0.02s
% Computational Cost: add. (57->14), mult. (16->4), div. (0->0), fcn. (40->4), ass. (0->10)
t31 = qJ(3) + qJ(4) + qJ(5);
t29 = sin(t31);
t32 = sin(qJ(1));
t35 = t32 * t29;
t30 = cos(t31);
t33 = cos(qJ(1));
t34 = t33 * t30;
t28 = t33 * t29;
t27 = t32 * t30;
t1 = [t28, 0, t27, t27, t27, 0; t35, 0, -t34, -t34, -t34, 0; 0, 0, -t29, -t29, -t29, 0; t34, 0, -t35, -t35, -t35, 0; t27, 0, t28, t28, t28, 0; 0, 0, -t30, -t30, -t30, 0; -t32, 0, 0, 0, 0, 0; t33, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
JR_rot  = t1;
